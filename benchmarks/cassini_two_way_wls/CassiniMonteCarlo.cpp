#include "CassiniMonteCarlo.hpp"
#include "simulation/CassiniRadioKernels.hpp"
#include <Eigen/Dense>
#include <array>
#include <cerrno>
#include <cmath>
#include <fcntl.h>
#include <iostream>
#include <regex>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char **environ;

namespace cassini_benchmark {
using namespace fd::simulation::monte_carlo;
namespace {
std::string number(double value) {
  std::ostringstream out;
  out << std::setprecision(17) << value;
  return out.str();
}
std::string read(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("Cannot read " + path.string());
  std::ostringstream out;
  out << in.rdbuf();
  return out.str();
}
std::string fingerprint(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("Cannot fingerprint " + path.string());
  std::uint64_t hash = 14695981039346656037ULL;
  char buffer[65536];
  while (in.read(buffer, sizeof buffer) || in.gcount())
    for (std::streamsize i = 0; i < in.gcount(); ++i) {
      hash ^= static_cast<unsigned char>(buffer[i]);
      hash *= 1099511628211ULL;
    }
  if (in.bad())
    throw std::runtime_error("Input fingerprint failed");
  std::ostringstream out;
  out << std::hex << hash;
  return out.str();
}
// No shell interpolation. Each fit gets a new process and independent SPICE
// state.
int execute(const std::vector<std::string> &arguments,
            const std::filesystem::path &log) {
  std::vector<char *> argv;
  for (const auto &value : arguments)
    argv.push_back(const_cast<char *>(value.c_str()));
  argv.push_back(nullptr);
  posix_spawn_file_actions_t actions;
  int code = posix_spawn_file_actions_init(&actions);
  if (code)
    throw std::runtime_error("Cannot initialize child process");
  code = posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, log.c_str(),
                                          O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (!code)
    code = posix_spawn_file_actions_adddup2(&actions, STDOUT_FILENO,
                                            STDERR_FILENO);
  pid_t pid = 0;
  if (!code)
    code = posix_spawn(&pid, argv[0], &actions, nullptr, argv.data(), environ);
  posix_spawn_file_actions_destroy(&actions);
  if (code)
    throw std::runtime_error("Cannot launch Cassini executable: " +
                             std::to_string(code));
  int status;
  while (waitpid(pid, &status, 0) < 0)
    if (errno != EINTR)
      throw std::runtime_error("Cannot wait for Cassini executable");
  return WIFEXITED(status) ? WEXITSTATUS(status)
                           : 128 + (WIFSIGNALED(status) ? WTERMSIG(status) : 0);
}
double scalar(const std::string &json, const std::string &key) {
  const std::regex pattern("\"" + key + "\"\\s*:\\s*([-+0-9.eE]+)");
  std::smatch match;
  if (!std::regex_search(json, match, pattern))
    throw std::runtime_error("Missing numeric summary field " + key);
  return std::stod(match[1]);
}
double nees(const Eigen::VectorXd &error, const Eigen::MatrixXd &covariance) {
  const Eigen::VectorXd sigma = covariance.diagonal().cwiseSqrt();
  if (!sigma.allFinite() || (sigma.array() <= 0).any())
    throw std::runtime_error("Invalid covariance scale");
  const Eigen::MatrixXd normalized = sigma.cwiseInverse().asDiagonal() *
                                     covariance *
                                     sigma.cwiseInverse().asDiagonal();
  const Eigen::LLT<Eigen::MatrixXd> factor(normalized);
  if (factor.info() != Eigen::Success)
    throw std::runtime_error("Covariance is not positive definite");
  const Eigen::VectorXd residual = error.cwiseQuotient(sigma);
  return residual.dot(factor.solve(residual));
}
} // namespace
void CassiniMonteCarloConfig::validate() const {
  if (!std::filesystem::is_regular_file(executable) || sourceRoot.empty() ||
      studyDirectory.empty() ||
      !(days >= 1 && days <= 14 && days == std::floor(days)) || vlbiCount < 1 ||
      vlbiCount > 1000 || maxIterations < 2 || maxIterations > 30)
    throw std::invalid_argument(
        "Invalid Cassini Monte Carlo paths/arc/iteration settings");
  for (double value : {rangeSigmaM, rangeBiasSigmaM, dopplerSigmaHz, vlbiSigmaM,
                       initialGuessScale})
    if (!std::isfinite(value) || value <= 0)
      throw std::invalid_argument("Invalid Cassini noise/guess scale");
}
CassiniWlsExperiment::CassiniWlsExperiment(CassiniMonteCarloConfig config)
    : config_(std::move(config)) {
  config_.validate();
}
std::vector<std::string>
CassiniWlsExperiment::arguments(const std::filesystem::path &output,
                                int count) const {
  return {config_.executable.string(),
          "--pass-schedule",
          "--days",
          number(config_.days),
          "--range-sigma-m",
          number(config_.rangeSigmaM),
          "--range-bias-sigma-m",
          number(config_.rangeBiasSigmaM),
          "--doppler-sigma-hz",
          number(config_.dopplerSigmaHz),
          "--vlbi-sigma-m",
          number(config_.vlbiSigmaM),
          "--initial-guess-scale",
          number(config_.initialGuessScale),
          "--max-iterations",
          std::to_string(config_.maxIterations),
          "--vlbi-count",
          std::to_string(count),
          "--output",
          output.string()};
}
std::string CassiniWlsExperiment::identity() const {
  std::ostringstream out;
  out << "Cassini fixed-truth experiment v1\n";
  for (const auto &option : arguments("OUTPUT", config_.vlbiCount))
    out << std::quoted(option) << ' ';
  out << '\n';
  std::vector<std::filesystem::path> inputs =
      fd::simulation::cassiniRadioKernels(config_.sourceRoot / "Kernels");
  for (const char *file :
       {"fk/cas_v40.tf", "sclk/cas00172.tsc", "ck/04256_04292py_as_flown.bc",
        "ck/04292_04320py_as_flown.bc", "ck/04281_04286ra.bc",
        "ck/04286_04291ra.bc", "ck/04291_04296ra.bc", "ck/04296_04301ra.bc",
        "ck/04293_04298ra.bc", "spk/200128RU_SCPSE_04251_04303.bsp.lbl"})
    inputs.push_back(config_.sourceRoot / "Kernels/Cassini" / file);
  for (const char *file : {"trajectory.csv", "metadata.json"})
    inputs.push_back(config_.sourceRoot / "results/cassini_two_way_wls/pass_schedule/truth" / file);
  inputs.push_back(config_.executable);
  for (const char *directory :
       {"include", "src", "benchmarks/cassini_two_way_wls",
        "benchmarks/cassini_trajectory"})
    for (const auto &entry : std::filesystem::recursive_directory_iterator(
             config_.sourceRoot / directory))
      if (entry.is_regular_file() && (entry.path().extension() == ".hpp" ||
                                      entry.path().extension() == ".cpp"))
        inputs.push_back(entry.path());
  std::sort(inputs.begin(), inputs.end());
  for (const auto &path : inputs)
    out << std::quoted(path.string()) << ' ' << fingerprint(path) << '\n';
  // FNV-1a is an accidental-change fingerprint, not a cryptographic digest.
  return out.str();
}
void CassiniWlsExperiment::prepare() {
  const auto cache = config_.studyDirectory / "prepared";
  std::filesystem::create_directories(cache);
  const auto cacheIdentity = [&] {
    std::ostringstream text;
    for (const char *file : {"observations.csv", "passes.csv",
                             "vlbi_observations.csv", "schedule.csv"})
      text << file << ' ' << fingerprint(cache / file) << '\n';
    return text.str();
  };
  if (!std::filesystem::exists(cache / "complete")) {
    auto args = arguments(cache, config_.vlbiCount);
    args.push_back("--generate-only");
    args.push_back("--no-noise");
    std::cout
        << "Preparing fixed Cassini geometry and noiseless observations.\n"
        << std::flush;
    if (execute(args, cache / "prepare.log") != 0)
      throw std::runtime_error("Cassini preparation failed; see " +
                               (cache / "prepare.log").string());
    {
      std::ofstream marker(cache / "complete.partial");
      marker << cacheIdentity();
      marker.close();
      if (!marker)
        throw std::runtime_error("Cannot save preparation marker");
    }
    std::filesystem::rename(cache / "complete.partial", cache / "complete");
  }
  if (read(cache / "complete") != cacheIdentity())
    throw std::runtime_error("Prepared observation cache changed");
  radio_ = std::make_unique<CsvTable>(cache / "observations.csv");
  passes_ = std::make_unique<CsvTable>(cache / "passes.csv");
  vlbi_ = std::make_unique<CsvTable>(cache / "vlbi_observations.csv");
}
CassiniTrialData
CassiniWlsExperiment::generateTrial(const TrialContext &context) {
  const auto directory = context.directory / "observations";
  std::filesystem::create_directories(directory);
  auto radio = *radio_, passes = *passes_, vlbi = *vlbi_;
  auto rangeRng = context.randomStream("range.white"),
       biasRng = context.randomStream("range.pass_bias"),
       dopplerRng = context.randomStream("doppler"),
       vlbiRng = context.randomStream("vlbi");
  std::normal_distribution<double> rangeNormal, biasNormal, dopplerNormal,
      vlbiNormal;
  std::map<int, double> biases;
  for (auto &row : passes.rows) {
    const int id = static_cast<int>(passes.number(row, "pass_id"));
    const double bias = config_.rangeBiasSigmaM * biasNormal(biasRng);
    biases[id] = bias;
    passes.set(row, "true_range_bias_m", bias);
  }
  for (auto &row : radio.rows) {
    const double doppler = config_.dopplerSigmaHz * dopplerNormal(dopplerRng);
    radio.set(row, "doppler_noise_hz", doppler);
    radio.set(row, "observed_doppler_hz",
              radio.number(row, "truth_doppler_hz") + doppler);
    const double bias =
        biases.at(static_cast<int>(radio.number(row, "pass_id")));
    radio.set(row, "range_bias_m", bias);
    if (radio.number(row, "range_active") == 1) {
      const double independent = config_.rangeSigmaM * rangeNormal(rangeRng);
      radio.set(row, "range_independent_noise_m", independent);
      radio.set(row, "range_noise_m", bias + independent);
      radio.set(row, "observed_range_km",
                radio.number(row, "truth_range_km") +
                    (bias + independent) / 1000);
    }
  }
  for (auto &row : vlbi.rows) {
    const double noise = config_.vlbiSigmaM * vlbiNormal(vlbiRng);
    vlbi.set(row, "noise_m", noise);
    vlbi.set(row, "observed_delay_km",
             vlbi.number(row, "truth_delay_km") + noise / 1000);
  }
  radio.write(directory / "observations.csv");
  passes.write(directory / "passes.csv");
  vlbi.write(directory / "vlbi_observations.csv");
  std::filesystem::copy_file(config_.studyDirectory / "prepared/schedule.csv",
                             directory / "schedule.csv",
                             std::filesystem::copy_options::overwrite_existing);
  std::ofstream seeds(directory / "streams.txt");
  for (const char *name : {"range.white", "range.pass_bias", "doppler", "vlbi"})
    seeds << name << ' ' << context.streamSeed(name) << '\n';
  return {directory};
}
TrialResult CassiniWlsExperiment::runVariant(const TrialContext &context,
                                             const CassiniTrialData &data,
                                             const std::string &variant) {
  if (variant != "no_vlbi" && variant != "with_vlbi")
    throw std::invalid_argument("Unknown Cassini estimator variant");
  const auto output = context.directory / variant;
  std::filesystem::create_directories(output);
  auto args = arguments(output, variant == "with_vlbi" ? config_.vlbiCount : 0);
  args.push_back("--observation-input");
  args.push_back(data.observations.string());
  std::cout << "Trial " << context.trialId << " / " << variant << "\n"
            << std::flush;
  const int code = execute(args, output / "fit.log");
  TrialResult result;
  if (code != 0 && code != 2) {
    result.status = TrialStatus::Failed;
    result.message =
        "Cassini fit exited " + std::to_string(code) + "; see fit.log";
    result.artifacts = {"fit.log"};
    return result;
  }
  result.status = code == 0 ? TrialStatus::Success : TrialStatus::NotConverged;
  result.message = code == 0 ? "" : "WLS did not converge";
  const auto summary = read(output / "summary.json");
  for (const auto &[key, unit] :
       std::vector<std::pair<std::string, std::string>>{
           {"fitted_position_rms_m", "m"},
           {"fitted_endpoint_error_m", "m"},
           {"postfit_range_rms_m", "m"},
           {"postfit_doppler_rms_hz", "Hz"},
           {"postfit_vlbi_rms_m", "m"},
           {"iterations", "count"},
           {"weighted_postfit_rms", "1"}})
    result.metrics[key] = {scalar(summary, key), unit};
  const CsvTable parameters(output / "initial_estimate.csv");
  for (const auto &row : parameters.rows)
    result.parameters.push_back(
        {row[parameters.column("component")], row[parameters.column("unit")],
         parameters.number(row, "truth"), parameters.number(row, "estimate")});
  Eigen::MatrixXd covariance(8, 8);
  std::ifstream cov(output / "initial_covariance.csv");
  std::string line, field;
  for (int i = 0; i < 8; ++i) {
    if (!std::getline(cov, line))
      throw std::runtime_error("Incomplete fit covariance");
    std::istringstream row(line);
    for (int j = 0; j < 8; ++j) {
      if (!std::getline(row, field, ','))
        throw std::runtime_error("Incomplete covariance row");
      const double value = std::stod(field);
      covariance(i, j) = value;
      result.covariance.push_back(value);
    }
  }
  Eigen::VectorXd error(8);
  for (int i = 0; i < 8; ++i)
    error[i] = result.parameters[i].estimate - result.parameters[i].truth;
  result.metrics["initial_parameter_nees"] = {nees(error, covariance), "1"};
  const CsvTable epochs(output / "epoch_error_covariance.csv");
  double total = 0;
  for (const auto &row : epochs.rows) {
    Eigen::VectorXd e(6);
    Eigen::MatrixXd p(6, 6);
    const std::array<const char *, 6> columns{
        "dx_m", "dy_m", "dz_m", "dvx_mm_s", "dvy_mm_s", "dvz_mm_s"};
    for (int i = 0; i < 6; ++i) {
      e[i] = epochs.number(row, columns[i]);
      for (int j = 0; j < 6; ++j)
        p(i, j) =
            epochs.number(row, "cov_" + std::to_string(i) + std::to_string(j));
    }
    total += nees(e, p);
  }
  result.metrics["mean_epoch_state_nees"] = {total / epochs.rows.size(), "1"};
  result.artifacts = {"summary.json", "initial_estimate.csv",
                      "initial_covariance.csv", "epoch_error_covariance.csv",
                      "fit.log"};
  const auto other = context.directory /
                     (variant == "with_vlbi" ? "no_vlbi" : "with_vlbi") /
                     "observations.csv";
  if (std::filesystem::exists(other) &&
      read(other) != read(output / "observations.csv"))
    throw std::runtime_error("Paired radio observations differ");
  if (context.retainDiagnostics)
    for (const char *file :
         {"trajectory.csv", "observations.csv", "residuals.csv",
          "vlbi_observations.csv", "vlbi_residuals.csv", "passes.csv",
          "schedule.csv", "iterations.csv"})
      result.artifacts.emplace_back(file);
  else
    for (const char *file :
         {"trajectory.csv", "residuals.csv", "vlbi_observations.csv",
          "vlbi_residuals.csv", "passes.csv", "schedule.csv", "iterations.csv"})
      std::filesystem::remove(output / file);
  // Always retain the small radio-input export so paired-data equality is
  // checked even after a metrics-only fit or a partial resume.
  if (!context.retainDiagnostics)
    result.artifacts.emplace_back("observations.csv");
  return result;
}
} // namespace cassini_benchmark
