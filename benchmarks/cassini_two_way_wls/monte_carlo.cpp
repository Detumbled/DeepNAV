#include "CassiniMonteCarlo.hpp"
#include <iostream>

int main(int argc, char **argv) {
  try {
    namespace mc = fd::simulation::monte_carlo;
    mc::MonteCarloConfig run;
    run.variants = {"no_vlbi", "with_vlbi"};
    run.outputDirectory = "results/cassini_two_way_wls/monte_carlo";
    cassini_benchmark::CassiniMonteCarloConfig experiment;
    experiment.sourceRoot = DEEPNAV_SOURCE_DIR;
    experiment.executable = std::filesystem::absolute(argv[0]).parent_path() /
                            "cassini_two_way_wls";
    for (int i = 1; i < argc; ++i) {
      const std::string key = argv[i];
      if (key == "--resume") {
        run.resume = true;
        continue;
      }
      if (key == "--stop-on-failure") {
        run.failurePolicy = mc::FailurePolicy::Stop;
        continue;
      }
      if (key == "--help") {
        std::cout
            << "--trials N --master-seed N --output DIR --resume "
               "--stop-on-failure --detail metrics|selected|full --days D "
               "--vlbi-count N --range-sigma-m M --range-bias-sigma-m M "
               "--doppler-sigma-hz H --vlbi-sigma-m M --initial-guess-scale F "
               "--max-iterations N --executable PATH\n";
        return 0;
      }
      if (++i == argc)
        throw std::invalid_argument("Missing value for " + key);
      const std::string value = argv[i];
      const auto integer = [&]() {
        if (value.empty() ||
            !std::all_of(value.begin(), value.end(),
                         [](char c) { return c >= '0' && c <= '9'; }))
          throw std::invalid_argument("Expected unsigned integer");
        return std::stoull(value);
      };
      if (key == "--trials")
        run.trialCount = integer();
      else if (key == "--master-seed")
        run.masterSeed = integer();
      else if (key == "--output")
        run.outputDirectory = value;
      else if (key == "--executable")
        experiment.executable = std::filesystem::absolute(value);
      else if (key == "--detail") {
        if (value == "metrics")
          run.detailLevel = mc::DetailLevel::Metrics;
        else if (value == "selected")
          run.detailLevel = mc::DetailLevel::Selected;
        else if (value == "full")
          run.detailLevel = mc::DetailLevel::Full;
        else
          throw std::invalid_argument("Invalid detail level");
      } else if (key == "--days")
        experiment.days = std::stod(value);
      else if (key == "--range-sigma-m")
        experiment.rangeSigmaM = std::stod(value);
      else if (key == "--range-bias-sigma-m")
        experiment.rangeBiasSigmaM = std::stod(value);
      else if (key == "--doppler-sigma-hz")
        experiment.dopplerSigmaHz = std::stod(value);
      else if (key == "--vlbi-sigma-m")
        experiment.vlbiSigmaM = std::stod(value);
      else if (key == "--initial-guess-scale")
        experiment.initialGuessScale = std::stod(value);
      else if (key == "--vlbi-count") {
        const auto n = integer();
        if (n > 1000)
          throw std::invalid_argument("VLBI count too large");
        experiment.vlbiCount = n;
      } else if (key == "--max-iterations") {
        const auto n = integer();
        if (n > 30)
          throw std::invalid_argument("Iteration limit too large");
        experiment.maxIterations = n;
      } else
        throw std::invalid_argument("Unknown option " + key);
    }
    run.outputDirectory = std::filesystem::absolute(run.outputDirectory);
    experiment.studyDirectory = run.outputDirectory;
    cassini_benchmark::CassiniWlsExperiment adapter(experiment);
    const auto result = mc::MonteCarloRunner(run).run(adapter);
    std::cout << "Executed " << result.executed << ", resumed "
              << result.resumed << ", successful " << result.succeeded
              << ", not converged " << result.notConverged << ", failed "
              << result.failed << '\n';
    return result.failed || result.notConverged ? 2 : 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
