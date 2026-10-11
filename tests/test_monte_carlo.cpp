#include "simulation/monte_carlo/MonteCarlo.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
using namespace fd::simulation::monte_carlo;
namespace {
void check(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}
template <class Function> void rejects(Function action, const char *message) {
  bool rejected = false;
  try {
    action();
  } catch (const std::exception &) {
    rejected = true;
  }
  check(rejected, message);
}
struct Data {
  double observation;
};
class Experiment : public MonteCarloExperiment<Data> {
public:
  int prepared = 0, generated = 0, fitted = 0;
  bool fail = false;
  std::string version = "linear experiment v1";
  std::string identity() const override { return version; }
  void prepare() override { ++prepared; }
  Data generateTrial(const TrialContext &context) override {
    ++generated;
    if (fail && context.trialId == 1)
      throw std::runtime_error("generation failure");
    auto rng = context.randomStream("sensor.white");
    return {std::normal_distribution<double>{}(rng)};
  }
  TrialResult runVariant(const TrialContext &, const Data &data,
                         const std::string &variant) override {
    ++fitted;
    const double scale = variant == "a" ? 1 : 2;
    TrialResult r;
    r.metrics.emplace("error", Metric{data.observation * scale, "m"});
    r.parameters.push_back({"x", "m", 0, data.observation * scale});
    r.covariance = {scale * scale};
    return r;
  }
};
} // namespace
int main() {
  std::filesystem::path root;
  try {
    root = std::filesystem::temp_directory_path() /
           ("deepnav-mc-test-" +
            std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()));
    MonteCarloConfig config;
    config.outputDirectory = root / "study";
    config.variants = {"a", "b"};
    Experiment experiment;
    const auto first = MonteCarloRunner(config).run(experiment);
    check(first.executed == 6 && first.succeeded == 6 &&
              experiment.prepared == 1 && experiment.generated == 3,
          "Preparation/trial counts incorrect");
    ResultWriter writer(config.outputDirectory);
    for (int id = 0; id < 3; ++id)
      check(writer.read(id, "b").metrics.at("error").value ==
                2 * writer.read(id, "a").metrics.at("error").value,
            "Paired data changed between variants");
    const TrialContext context(2026, 0, root, false),
        other(2026, 1, root, false);
    auto rng = context.randomStream("sensor.white");
    const auto value = rng();
    context.randomStream("new.sensor");
    auto replay = context.randomStream("sensor.white");
    check(value == replay() &&
              context.streamSeed("sensor.white") !=
                  other.streamSeed("sensor.white") &&
              context.streamSeed("sensor.white") !=
                  context.streamSeed("new.sensor"),
          "Stream isolation/replay failed");
    config.resume = true;
    Experiment resumed;
    const auto second = MonteCarloRunner(config).run(resumed);
    check(second.resumed == 6 && second.executed == 0 && resumed.prepared == 0,
          "Resume reran completed work");
    // Simulate interruption before one variant's atomic commit.
    std::filesystem::remove(writer.variantDirectory(1, "b") / "result.txt");
    std::ofstream(writer.variantDirectory(1, "b") / "result.txt.partial")
        << "interrupted";
    const auto partial = MonteCarloRunner(config).run(resumed);
    check(partial.executed == 1 && partial.resumed == 5 &&
              resumed.generated == 1,
          "Partial resume did not isolate pending variant");
    check(writer.read(1, "b").metrics.at("error").value ==
              2 * writer.read(1, "a").metrics.at("error").value,
          "Partial replay changed observations");
    Experiment changed;
    changed.version = "changed inputs";
    rejects([&] { MonteCarloRunner(config).run(changed); },
            "Changed identity accepted on resume");
    auto different = config;
    different.masterSeed++;
    rejects([&] { MonteCarloRunner(different).run(resumed); },
            "Changed seed accepted on resume");
    config.outputDirectory = root / "repeat";
    config.resume = false;
    Experiment repeat;
    MonteCarloRunner(config).run(repeat);
    ResultWriter repeated(config.outputDirectory);
    check(repeated.read(2, "a").metrics.at("error").value ==
              writer.read(2, "a").metrics.at("error").value,
          "Full replay differs");
    rejects([&] { MonteCarloRunner(config).run(repeat); },
            "Existing study overwritten");
    config.outputDirectory = root / "failures";
    Experiment failing;
    failing.fail = true;
    const auto failures = MonteCarloRunner(config).run(failing);
    check(failures.failed == 2 && failures.succeeded == 4 &&
              failures.executed == 6,
          "Failed generation disappeared from denominator");
    config.outputDirectory = root / "stop";
    config.failurePolicy = FailurePolicy::Stop;
    rejects([&] { MonteCarloRunner(config).run(failing); },
            "Stop policy ignored");
    check(std::filesystem::exists(config.outputDirectory / "metrics.csv") &&
              !std::filesystem::exists(config.outputDirectory / "trial_2"),
          "Stop did not persist outcomes");
    auto invalid = config;
    invalid.variants = {"../escape"};
    rejects([&] { invalid.validate(); }, "Unsafe variant path accepted");
    TrialResult bad;
    bad.metrics["nan"] = {std::numeric_limits<double>::quiet_NaN(), "m"};
    rejects([&] { bad.validate(); }, "Nonfinite metric accepted");
    std::ofstream(writer.variantDirectory(2, "a") / "result.txt")
        << "MC_RESULT_1 0 \"\"\n1\n";
    rejects([&] { writer.read(2, "a"); },
            "Truncated committed result accepted");
    std::filesystem::remove_all(root);
    std::cout << "Monte Carlo replay, pairing, resume and failure policies "
                 "verified.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    if (!root.empty())
      std::filesystem::remove_all(root);
    return 1;
  }
}
