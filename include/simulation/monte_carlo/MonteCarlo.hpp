#pragma once

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <map>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace fd::simulation::monte_carlo {

enum class FailurePolicy { Continue, Stop };
enum class DetailLevel { Metrics, Selected, Full };
enum class TrialStatus { Success, NotConverged, Failed };

struct MonteCarloConfig {
  std::size_t trialCount = 3;
  std::uint64_t masterSeed = 2026;
  std::vector<std::string> variants;
  std::filesystem::path outputDirectory;
  bool resume = false;
  FailurePolicy failurePolicy = FailurePolicy::Continue;
  DetailLevel detailLevel = DetailLevel::Selected;
  std::vector<std::size_t> diagnosticTrialIds{0};
  void validate() const;
};

// Stateless seed derivation: requesting another stream never consumes this one.
// Engine sequences are portable; std::normal_distribution transforms need not
// be identical between standard-library implementations.
class TrialContext {
  std::uint64_t masterSeed_;

public:
  const std::size_t trialId;
  const std::filesystem::path directory;
  const bool retainDiagnostics;
  TrialContext(std::uint64_t seed, std::size_t id, std::filesystem::path path,
               bool diagnostics);
  std::uint64_t streamSeed(const std::string &name) const;
  std::mt19937_64 randomStream(const std::string &name) const;
};

struct Metric {
  double value = 0;
  std::string unit;
};
struct ParameterEstimate {
  std::string name, unit;
  double truth = 0, estimate = 0;
};
struct TrialResult {
  TrialStatus status = TrialStatus::Success;
  std::string message;
  std::map<std::string, Metric> metrics;
  std::vector<ParameterEstimate> parameters;
  // Row-major covariance in the same component order and units as parameters.
  std::vector<double> covariance;
  // Paths relative to the variant directory. Artifacts are experiment-owned.
  std::vector<std::filesystem::path> artifacts;
  void validate() const;
};

class ResultWriter {
  std::filesystem::path root_;

public:
  explicit ResultWriter(std::filesystem::path root);
  void initialize(const MonteCarloConfig &,
                  const std::string &experimentIdentity);
  std::filesystem::path variantDirectory(std::size_t trial,
                                         const std::string &variant) const;
  bool completed(std::size_t trial, const std::string &variant) const;
  TrialResult read(std::size_t trial, const std::string &variant) const;
  void write(std::size_t trial, const std::string &variant,
             const TrialResult &);
  void aggregate(const MonteCarloConfig &);
};

// Typed trial payloads keep experiment-specific observations and truth out of
// the generic runner. Variants receive the same immutable payload.
template <class TrialData> class MonteCarloExperiment {
public:
  virtual ~MonteCarloExperiment() = default;
  // Include canonical experiment options and input/source fingerprints.
  virtual std::string identity() const = 0;
  virtual void prepare() = 0;
  virtual TrialData generateTrial(const TrialContext &) = 0;
  virtual TrialResult runVariant(const TrialContext &, const TrialData &,
                                 const std::string &variant) = 0;
};

struct RunSummary {
  std::size_t executed = 0, resumed = 0, succeeded = 0, notConverged = 0,
              failed = 0;
};

class MonteCarloRunner {
  MonteCarloConfig config_;

public:
  explicit MonteCarloRunner(MonteCarloConfig config);
  template <class TrialData>
  RunSummary run(MonteCarloExperiment<TrialData> &experiment) const {
    ResultWriter writer(config_.outputDirectory);
    writer.initialize(config_, experiment.identity());
    RunSummary summary;
    bool prepared = false;
    for (std::size_t id = 0; id < config_.trialCount; ++id) {
      bool pending = false;
      for (const auto &variant : config_.variants)
        pending |= !writer.completed(id, variant);
      const bool diagnostics =
          config_.detailLevel == DetailLevel::Full ||
          (config_.detailLevel == DetailLevel::Selected &&
           std::find(config_.diagnosticTrialIds.begin(),
                     config_.diagnosticTrialIds.end(),
                     id) != config_.diagnosticTrialIds.end());
      const TrialContext context(config_.masterSeed, id,
                                 config_.outputDirectory /
                                     ("trial_" + std::to_string(id)),
                                 diagnostics);
      const auto record = [&](const TrialResult &result) {
        switch (result.status) {
        case TrialStatus::Success:
          ++summary.succeeded;
          break;
        case TrialStatus::NotConverged:
          ++summary.notConverged;
          break;
        case TrialStatus::Failed:
          ++summary.failed;
          break;
        }
      };
      const auto runPending = [&](const TrialData *data,
                                  const std::string &generationFailure) {
        for (const auto &variant : config_.variants) {
          if (writer.completed(id, variant)) {
            ++summary.resumed;
            record(writer.read(id, variant));
            continue;
          }
          TrialResult result;
          try {
            if (!data)
              throw std::runtime_error(generationFailure);
            result = experiment.runVariant(context, *data, variant);
            result.validate();
          } catch (const std::exception &error) {
            result = TrialResult{};
            result.status = TrialStatus::Failed;
            result.message = error.what();
          }
          writer.write(id, variant, result);
          ++summary.executed;
          record(result);
          if (result.status != TrialStatus::Success &&
              config_.failurePolicy == FailurePolicy::Stop) {
            writer.aggregate(config_);
            throw std::runtime_error("Monte Carlo stopped at trial " +
                                     std::to_string(id) + ", " + variant +
                                     ": " + result.message);
          }
        }
      };
      if (!pending) {
        runPending(nullptr, "");
        continue;
      }
      if (!prepared) {
        experiment.prepare();
        prepared = true;
      }
      std::optional<TrialData> data;
      std::string failure;
      try {
        data.emplace(experiment.generateTrial(context));
      } catch (const std::exception &error) {
        failure = error.what();
      }
      runPending(data ? &*data : nullptr, failure);
    }
    writer.aggregate(config_);
    return summary;
  }
};
} // namespace fd::simulation::monte_carlo
