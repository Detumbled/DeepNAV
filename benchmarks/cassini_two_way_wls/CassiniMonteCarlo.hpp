#pragma once
#include "CsvTable.hpp"
#include "simulation/monte_carlo/MonteCarlo.hpp"
#include <memory>

namespace cassini_benchmark {
struct CassiniMonteCarloConfig {
  std::filesystem::path executable, sourceRoot, studyDirectory;
  double days = 14, rangeSigmaM = 1, rangeBiasSigmaM = 10, dopplerSigmaHz = .01;
  double vlbiSigmaM = .01, initialGuessScale = 10;
  int vlbiCount = 60, maxIterations = 10;
  void validate() const;
};
struct CassiniTrialData {
  std::filesystem::path observations;
};
class CassiniWlsExperiment final
    : public fd::simulation::monte_carlo::MonteCarloExperiment<
          CassiniTrialData> {
  CassiniMonteCarloConfig config_;
  std::unique_ptr<CsvTable> radio_, passes_, vlbi_;
  std::vector<std::string> arguments(const std::filesystem::path &output,
                                     int vlbiCount) const;

public:
  explicit CassiniWlsExperiment(CassiniMonteCarloConfig);
  std::string identity() const override;
  void prepare() override;
  CassiniTrialData
  generateTrial(const fd::simulation::monte_carlo::TrialContext &) override;
  fd::simulation::monte_carlo::TrialResult
  runVariant(const fd::simulation::monte_carlo::TrialContext &,
             const CassiniTrialData &, const std::string &) override;
};
} // namespace cassini_benchmark
