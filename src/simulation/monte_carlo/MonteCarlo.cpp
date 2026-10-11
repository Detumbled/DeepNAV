#include "simulation/monte_carlo/MonteCarlo.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>

namespace fd::simulation::monte_carlo {
namespace {
bool safeName(const std::string &name) {
  return !name.empty() && name != "." && name != ".." &&
         std::all_of(name.begin(), name.end(), [](unsigned char c) {
           return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                  (c >= '0' && c <= '9') || c == '_' || c == '-';
         });
}
std::string csv(const std::string &value) {
  std::string out = "\"";
  for (char c : value) {
    if (c == '\"')
      out += '\"';
    out += c;
  }
  return out + '\"';
}
void atomicWrite(const std::filesystem::path &path,
                 const std::string &contents) {
  std::filesystem::create_directories(path.parent_path());
  const auto temporary = path.string() + ".partial";
  {
    std::ofstream out(temporary, std::ios::binary);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << contents;
    out.close();
  }
  std::filesystem::rename(temporary, path);
}
std::string contents(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input)
    throw std::runtime_error("Cannot read " + path.string());
  std::ostringstream out;
  out << input.rdbuf();
  return out.str();
}
std::uint64_t mix(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}
} // namespace
void MonteCarloConfig::validate() const {
  if ((failurePolicy != FailurePolicy::Continue &&
       failurePolicy != FailurePolicy::Stop) ||
      (detailLevel != DetailLevel::Metrics &&
       detailLevel != DetailLevel::Selected &&
       detailLevel != DetailLevel::Full))
    throw std::invalid_argument("Invalid Monte Carlo policy");
  if (trialCount == 0 || variants.empty() || outputDirectory.empty())
    throw std::invalid_argument(
        "Monte Carlo needs trials, variants and an output directory");
  std::set<std::string> names;
  for (const auto &name : variants)
    if (!safeName(name) || !names.insert(name).second)
      throw std::invalid_argument("Invalid or duplicate variant name");
  std::set<std::size_t> ids;
  for (auto id : diagnosticTrialIds)
    if (id >= trialCount || !ids.insert(id).second)
      throw std::invalid_argument("Invalid or duplicate diagnostic trial ID");
}
TrialContext::TrialContext(std::uint64_t seed, std::size_t id,
                           std::filesystem::path path, bool diagnostics)
    : masterSeed_(seed), trialId(id), directory(std::move(path)),
      retainDiagnostics(diagnostics) {}
std::uint64_t TrialContext::streamSeed(const std::string &name) const {
  if (name.empty())
    throw std::invalid_argument("Random stream name is empty");
  std::uint64_t hash = 14695981039346656037ULL;
  for (unsigned char c : name) {
    hash ^= c;
    hash *= 1099511628211ULL;
  }
  return mix(mix(masterSeed_) ^ mix(static_cast<std::uint64_t>(trialId)) ^
             mix(hash));
}
std::mt19937_64 TrialContext::randomStream(const std::string &name) const {
  return std::mt19937_64(streamSeed(name));
}
void TrialResult::validate() const {
  if (status != TrialStatus::Success && status != TrialStatus::NotConverged &&
      status != TrialStatus::Failed)
    throw std::invalid_argument("Invalid trial status");
  for (const auto &[name, metric] : metrics)
    if (name.empty() || !std::isfinite(metric.value))
      throw std::invalid_argument("Invalid metric");
  std::set<std::string> names;
  for (const auto &p : parameters)
    if (p.name.empty() || !names.insert(p.name).second ||
        !std::isfinite(p.truth) || !std::isfinite(p.estimate))
      throw std::invalid_argument("Invalid parameter estimate");
  const auto n = parameters.size();
  if (!covariance.empty()) {
    if (covariance.size() != n * n)
      throw std::invalid_argument(
          "Covariance dimensions do not match parameters");
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j) {
        const double a = covariance[i * n + j], b = covariance[j * n + i];
        if (!std::isfinite(a) ||
            std::abs(a - b) >
                1e-10 * std::max({1e-30, std::abs(a), std::abs(b)}) ||
            (i == j && a < 0))
          throw std::invalid_argument("Invalid covariance");
      }
  }
  for (const auto &path : artifacts) {
    if (path.empty() || path.is_absolute())
      throw std::invalid_argument("Artifact paths must be relative");
    for (const auto &part : path)
      if (part == "..")
        throw std::invalid_argument("Artifact path escapes its directory");
  }
}
ResultWriter::ResultWriter(std::filesystem::path root)
    : root_(std::move(root)) {}
void ResultWriter::initialize(const MonteCarloConfig &config,
                              const std::string &identity) {
  config.validate();
  if (identity.empty())
    throw std::invalid_argument("Experiment identity is empty");
  std::ostringstream canonical;
  canonical << "DeepNAV Monte Carlo schema 1\n"
            << config.trialCount << ' ' << config.masterSeed << ' '
            << static_cast<int>(config.detailLevel) << ' '
            << static_cast<int>(config.failurePolicy) << '\n';
  for (const auto &variant : config.variants)
    canonical << std::quoted(variant) << ' ';
  canonical << '\n';
  for (auto id : config.diagnosticTrialIds)
    canonical << id << ' ';
  canonical << '\n' << std::quoted(identity) << '\n';
  const auto manifest = root_ / "manifest.txt";
  if (std::filesystem::exists(manifest)) {
    if (!config.resume)
      throw std::runtime_error("Study already exists; explicitly enable resume "
                               "or choose a new output directory");
    if (contents(manifest) != canonical.str())
      throw std::runtime_error(
          "Resume rejected: study configuration or experiment inputs changed");
  } else {
    if (std::filesystem::exists(root_) && !std::filesystem::is_empty(root_))
      throw std::runtime_error(
          "Output directory contains an unrecognized study");
    atomicWrite(manifest, canonical.str());
  }
}
std::filesystem::path
ResultWriter::variantDirectory(std::size_t trial,
                               const std::string &variant) const {
  if (!safeName(variant))
    throw std::invalid_argument("Invalid variant path");
  return root_ / ("trial_" + std::to_string(trial)) / variant;
}
bool ResultWriter::completed(std::size_t trial,
                             const std::string &variant) const {
  if (!std::filesystem::exists(variantDirectory(trial, variant) / "result.txt"))
    return false;
  read(trial, variant);
  return true;
}
void ResultWriter::write(std::size_t trial, const std::string &variant,
                         const TrialResult &result) {
  result.validate();
  std::ostringstream out;
  out << std::setprecision(17) << "MC_RESULT_1 "
      << static_cast<int>(result.status) << ' ' << std::quoted(result.message)
      << '\n';
  out << result.metrics.size() << '\n';
  for (const auto &[name, m] : result.metrics)
    out << std::quoted(name) << ' ' << m.value << ' ' << std::quoted(m.unit)
        << '\n';
  out << result.parameters.size() << '\n';
  for (const auto &p : result.parameters)
    out << std::quoted(p.name) << ' ' << std::quoted(p.unit) << ' ' << p.truth
        << ' ' << p.estimate << '\n';
  out << result.covariance.size() << '\n';
  for (double value : result.covariance)
    out << value << ' ';
  out << '\n';
  out << result.artifacts.size() << '\n';
  for (const auto &p : result.artifacts)
    out << std::quoted(p.generic_string()) << '\n';
  out << "END\n";
  atomicWrite(variantDirectory(trial, variant) / "result.txt", out.str());
}
TrialResult ResultWriter::read(std::size_t trial,
                               const std::string &variant) const {
  std::istringstream in(
      contents(variantDirectory(trial, variant) / "result.txt"));
  std::string tag;
  int status = 0;
  TrialResult result;
  in >> tag >> status >> std::quoted(result.message);
  if (tag != "MC_RESULT_1")
    throw std::runtime_error("Invalid result schema");
  result.status = static_cast<TrialStatus>(status);
  const auto count = [&] {
    std::size_t n = 0;
    in >> n;
    if (!in || n > 1000000)
      throw std::runtime_error("Invalid result length");
    return n;
  };
  auto n = count();
  for (std::size_t i = 0; i < n; ++i) {
    std::string name;
    Metric m;
    in >> std::quoted(name) >> m.value >> std::quoted(m.unit);
    if (!result.metrics.emplace(name, m).second)
      throw std::runtime_error("Duplicate stored metric");
  }
  n = count();
  for (std::size_t i = 0; i < n; ++i) {
    ParameterEstimate p;
    in >> std::quoted(p.name) >> std::quoted(p.unit) >> p.truth >> p.estimate;
    result.parameters.push_back(p);
  }
  n = count();
  result.covariance.resize(n);
  for (double &value : result.covariance)
    in >> value;
  n = count();
  for (std::size_t i = 0; i < n; ++i) {
    std::string path;
    in >> std::quoted(path);
    result.artifacts.emplace_back(path);
  }
  in >> tag;
  if (!in || tag != "END")
    throw std::runtime_error("Incomplete result file");
  in >> std::ws;
  if (!in.eof())
    throw std::runtime_error("Trailing result data");
  result.validate();
  return result;
}
void ResultWriter::aggregate(const MonteCarloConfig &config) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "trial_id,variant,status,metric,value,unit,message\n";
  for (std::size_t id = 0; id < config.trialCount; ++id)
    for (const auto &variant : config.variants) {
      if (!completed(id, variant))
        continue;
      const auto r = read(id, variant);
      const auto prefix =
          std::to_string(id) + "," + csv(variant) + "," +
          csv(r.status == TrialStatus::Success        ? "success"
              : r.status == TrialStatus::NotConverged ? "not_converged"
                                                      : "failed") +
          ",";
      if (r.metrics.empty())
        out << prefix << ",,," << csv(r.message) << '\n';
      for (const auto &[name, m] : r.metrics)
        out << prefix << csv(name) << ',' << m.value << ',' << csv(m.unit)
            << ',' << csv(r.message) << '\n';
    }
  atomicWrite(root_ / "metrics.csv", out.str());
}
MonteCarloRunner::MonteCarloRunner(MonteCarloConfig config)
    : config_(std::move(config)) {
  config_.validate();
}
} // namespace fd::simulation::monte_carlo
