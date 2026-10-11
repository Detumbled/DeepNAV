#include "../cassini_trajectory/SavedSimulation.hpp"
#include "../cassini_two_way_wls/CsvTable.hpp"
#include "Clocks/ClockTruthSimulator.hpp"
#include "Clocks/DSAC.hpp"
#include "ReceptionTrajectory.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "filters/CartesianPropagator.hpp"
#include "filters/EKF.hpp"
#include "observations/radiometric/OneWayMeasurement.hpp"
#include "observations/radiometric/PropagationCorrections.hpp"
#include "observations/radiometric/SpiceLightTime.hpp"
#include "simulation/CassiniRadioKernels.hpp"
#include <Eigen/Cholesky>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <random>

namespace radio = fd::observations::radiometric;
namespace clocks = fd::clocks;
namespace filters = fd::filters;
namespace fs = std::filesystem;
using Vector = Eigen::VectorXd;
struct Options {
  fs::path output = "results/cassini_one_way_ekf";
  std::string clock = "both";
  double days = 14, maxStep = 60, rangeSigmaM = 10, dopplerSigmaHz = .01,
         guessScale = 1, clockFrequencySigma = 1e-10;
  unsigned seed = 2026;
  bool noise = true;
};
Options parse(int argc, char **argv) {
  Options o;
  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    if (key == "--no-noise") {
      o.noise = false;
      continue;
    }
    if (key == "--help") {
      std::cout << "--output DIR --days D --clock both|ideal|dsac --max-step-s "
                   "S --range-sigma-m M --doppler-sigma-hz H "
                   "--initial-guess-scale F --clock-frequency-sigma Y --seed N "
                   "--no-noise\n";
      std::exit(0);
    }
    if (++i == argc)
      throw std::invalid_argument("Missing option value");
    if (key == "--output")
      o.output = argv[i];
    else if (key == "--clock")
      o.clock = argv[i];
    else {
      const double v = std::stod(argv[i]);
      if (!std::isfinite(v))
        throw std::invalid_argument("Nonfinite option");
      if (key == "--days")
        o.days = v;
      else if (key == "--max-step-s")
        o.maxStep = v;
      else if (key == "--range-sigma-m")
        o.rangeSigmaM = v;
      else if (key == "--doppler-sigma-hz")
        o.dopplerSigmaHz = v;
      else if (key == "--clock-frequency-sigma")
        o.clockFrequencySigma = v;
      else if (key == "--initial-guess-scale")
        o.guessScale = v;
      else if (key == "--seed" && v >= 0 && v <= 4294967295. &&
               v == std::floor(v))
        o.seed = static_cast<unsigned>(v);
      else
        throw std::invalid_argument("Unknown/invalid option " + key);
    }
  }
  if (o.days <= 0 || o.days > 14 || o.days != std::floor(o.days) ||
      o.maxStep <= 0 || o.maxStep > 60 || o.rangeSigmaM <= 0 ||
      o.dopplerSigmaHz <= 0 || o.guessScale <= 0 ||
      o.clockFrequencySigma <= 0 ||
      (o.clock != "both" && o.clock != "ideal" && o.clock != "dsac"))
    throw std::invalid_argument("Invalid configuration");
  return o;
}
struct ClockHistory {
  std::vector<clocks::ClockState> nodes;
  double epoch, step = 60, drift;
  clocks::ClockState stateAt(double t) const {
    const double u = (t - epoch) / step;
    if (u < 0 || u > nodes.size() - 1)
      throw std::out_of_range("Clock truth outside history");
    const auto i = std::min(static_cast<std::size_t>(u), nodes.size() - 2);
    const double dt = (u - i) * step;
    // Interpolate the Wiener bias increment; deterministic aging is exact.
    const double stochastic = nodes[i + 1].bias_s - nodes[i].bias_s -
                              nodes[i].fractional_frequency * step -
                              .5 * drift * step * step;
    return {nodes[i].bias_s + nodes[i].fractional_frequency * dt +
                .5 * drift * dt * dt + stochastic * dt / step,
            nodes[i].fractional_frequency + drift * dt};
  }
};
struct Pass {
  int id;
  std::size_t station;
  double begin, end;
};
struct Observation {
  double elapsed, elevation;
  std::size_t station;
  int pass;
  bool range;
  long double tag;
  Eigen::Vector2d truth, noise, value;
};
std::ofstream csv(const fs::path &path) {
  std::ofstream out(path);
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out << std::setprecision(17);
  return out;
}

int main(int argc, char **argv) {
  try {
    const auto options = parse(argc, argv);
    const fs::path root(DEEPNAV_SOURCE_DIR);
    od::SpiceErrorModeGuard errors("RETURN", "NONE");
    auto paths = fd::simulation::cassiniRadioKernels(root / "Kernels");
    for (const char *file :
         {"Cassini/fk/cas_v40.tf", "Cassini/sclk/cas00172.tsc",
          "Cassini/ck/04256_04292py_as_flown.bc",
          "Cassini/ck/04292_04320py_as_flown.bc", "Cassini/ck/04281_04286ra.bc",
          "Cassini/ck/04286_04291ra.bc", "Cassini/ck/04291_04296ra.bc",
          "Cassini/ck/04296_04301ra.bc", "Cassini/ck/04293_04298ra.bc"})
      paths.push_back(root / "Kernels" / file);
    fd::dynamics::SpiceKernelSet kernels(paths);
    const CassiniAttitude attitude(paths);
    const Sat427Gravity gravity(
        root / "Kernels/Cassini/spk/200128RU_SCPSE_04251_04303.bsp.lbl");
    double recoilEpoch;
    str2et_c("2005-08-13T00:00:00", &recoilEpoch);
    od::throwIfSpiceFailed("Thermal epoch");
    double radii[3];
    SpiceInt n;
    bodvrd_c("SATURN", "RADII", 3, &n, radii);
    od::throwIfSpiceFailed("Saturn radii");
    const fd::simulation::CassiniForceModel force(
        gravity, attitude, fd::simulation::cassiniPerturbers(gravity), 2, 3220,
        12, 1.42, {-2.725e-13, -2.725e-13, -5.232237415805135e-12}, recoilEpoch,
        radii[0]);
    fd::dynamics::SpiceStateProvider saturn("SATURN"), sun("SUN");
    const fd::simulation::SavedSimulation truth(
        root / "results/cassini_two_way_wls/pass_schedule/truth/trajectory.csv",
        force, saturn);
    const double t0 = truth.startEpoch(), duration = options.days * 86400;
    if (t0 + duration > truth.endEpoch())
      throw std::invalid_argument("Arc exceeds retained truth");
    const std::array<std::string, 3> names{"DSS-14", "DSS-43", "DSS-63"};
    const cassini_benchmark::CsvTable passInput(
        root /
        "results/cassini_two_way_wls/pass_schedule/with_vlbi/passes.csv");
    std::vector<Pass> passes;
    for (const auto &row : passInput.rows) {
      const auto name = row[passInput.column("station")];
      auto found = std::find(names.begin(), names.end(), name);
      if (found == names.end())
        throw std::invalid_argument("Unknown pass station");
      passes.push_back({static_cast<int>(passInput.number(row, "pass_id")),
                        static_cast<std::size_t>(found - names.begin()),
                        passInput.number(row, "start_elapsed_s"),
                        passInput.number(row, "end_elapsed_s")});
    }
    const double muSun = fd::simulation::cassiniKernelGM("SUN");
    const radio::SpiceTdbClock tdbMinusTt;
    const radio::ClockOffset nominalOffset = [&](double t) {
      return tdbMinusTt(t);
    };
    const radio::EndpointClock ground{t0, {}, nominalOffset};
    const radio::FrequencyRamp ramp{t0, 7.175e9, 0};
    std::vector<std::unique_ptr<fd::dynamics::SpiceStateProvider>> stations;
    std::vector<radio::OneWayLink> links;
    for (const auto &name : names) {
      stations.push_back(
          std::make_unique<fd::dynamics::SpiceStateProvider>(name));
      const auto delay = [&, name](const radio::LegSolution &leg, const auto &,
                                   const auto &) {
        const double elevation = radio::spiceStationElevationRadians(
            name, leg.emissionTdb,
            leg.receiver.positionKm - leg.emitter.positionKm);
        if (elevation <= 0)
          throw std::runtime_error("Uplink below horizon");
        return radio::pointMassShapiroSeconds(leg, sun, muSun) +
               radio::zenithPathDelaySeconds(2.3, elevation);
      };
      links.emplace_back(radio::libraryReceptionSolver({}, delay));
    }
    const fd::perturbations::AccelerationFunction acceleration =
        [&](double t, const auto &r, const auto &v) {
          return static_cast<fd::perturbations::AccelerationEvaluation>(
              force.evaluate(t, r, v, true));
        };
    for (const std::string kind : {"ideal", "dsac"}) {
      if (options.clock != "both" && options.clock != kind)
        continue;
      const auto wall = std::chrono::steady_clock::now();
      const bool ideal = kind == "ideal";
      const auto clock = ideal ? clocks::DSAC::fromParameters({})
                               : clocks::DSAC::shortTermWhiteFmBaseline();
      const auto parameters = clock.parameters();
      clocks::ClockTruthSimulator simulator(
          clock, static_cast<std::uint64_t>(options.seed) + 700001);
      ClockHistory clockTruth{{}, t0, 60, parameters.frequency_drift_per_s};
      clockTruth.nodes.push_back(ideal ? clocks::ClockState{}
                                       : clocks::ClockState{2e-6, 2e-12});
      for (int i = 0; i < static_cast<int>(duration / 60); ++i)
        clockTruth.nodes.push_back(simulator.step(clockTruth.nodes.back(), 60));
      const radio::EndpointClock receiverTruth{
          t0, [&](double t) { return clockTruth.stateAt(t).bias_s; },
          nominalOffset};
      std::mt19937_64 rangeRng(static_cast<std::uint64_t>(options.seed) + 1),
          dopplerRng(static_cast<std::uint64_t>(options.seed) + 2);
      std::normal_distribution<double> rangeGaussian, dopplerGaussian;
      std::vector<Observation> observations;
      int rangeCount = 0, rejected = 0;
      for (const auto &pass : passes) {
        for (double elapsed = pass.begin + 60;
             elapsed <= std::min(pass.end, duration - 60); elapsed += 60) {
          const auto j = pass.station;
          const double t = t0 + elapsed;
          // Geometric visibility first, so media never evaluates below the
          // horizon.
          const radio::OneWayLink geometry(radio::libraryReceptionSolver());
          const auto a = geometry.solve(t - 60, *stations[j], truth),
                     b = geometry.solve(t, *stations[j], truth);
          const auto elevation = [&](const auto &leg) {
            return radio::spiceStationElevationRadians(
                names[j], leg.emissionTdb,
                leg.receiver.positionKm - leg.emitter.positionKm);
          };
          const double e = std::min(elevation(a), elevation(b));
          if (e < 10 * std::numbers::pi / 180) {
            ++rejected;
            continue;
          }
          const long double tag = receiverTruth.reading(t);
          const auto counted =
              links[j].count(tag, 60, *stations[j], truth, ground,
                             receiverTruth, ramp, ramp.frequencyHz);
          const double actualElevation =
              std::min(elevation(counted.start), elevation(counted.end));
          if (actualElevation < 10 * std::numbers::pi / 180) {
            ++rejected;
            continue;
          }
          const bool range =
              static_cast<int>(std::llround(elapsed - pass.begin)) % 1200 == 0;
          Eigen::Vector2d noise{options.noise ? options.rangeSigmaM / 1000 *
                                                    rangeGaussian(rangeRng)
                                              : 0,
                                options.noise ? options.dopplerSigmaHz *
                                                    dopplerGaussian(dopplerRng)
                                              : 0};
          Eigen::Vector2d exact{counted.pseudorangeKm, counted.residualHz};
          observations.push_back(
              {elapsed, actualElevation * 180 / std::numbers::pi, j, pass.id,
               range, tag, exact, noise, exact + noise});
          if (range)
            ++rangeCount;
        }
      }
      if (observations.size() < 10 || !rangeCount)
        throw std::runtime_error("Insufficient one-way tracking");
      const fs::path folder = options.output / kind;
      fs::create_directories(folder);
      auto obsFile = csv(folder / "observations.csv");
      obsFile
          << "truth_elapsed_s,raw_receiver_tag_s,station,pass_id,range_active,"
             "minimum_elevation_deg,truth_range_km,truth_doppler_hz,range_"
             "noise_m,doppler_noise_hz,observed_range_km,observed_doppler_hz\n";
      for (const auto &a : observations)
        obsFile << a.elapsed << ',' << a.tag << ',' << names[a.station] << ','
                << a.pass << ',' << a.range << ',' << a.elevation << ','
                << a.truth[0] << ',' << a.truth[1] << ',' << a.noise[0] * 1000
                << ',' << a.noise[1] << ',' << a.value[0] << ',' << a.value[1]
                << '\n';
      filters::CartesianPropagationConfig config;
      config.clock = parameters;
      config.integrator.absoluteTolerance = 1e-13;
      config.integrator.relativeTolerance = 1e-13;
      config.integrator.initialStep = 10;
      config.integrator.maximumStep = options.maxStep;
      const filters::CartesianPropagator propagate(acceleration, config);
      filters::EKF filter(filters::StateLayout::OrbitClock);
      Vector initial(8);
      initial.head<6>() = truth.relative(t0);
      initial.tail<2>().setZero();
      initial.head<3>() += options.guessScale * Eigen::Vector3d(1, -.8, .5);
      initial.segment<3>(3) +=
          options.guessScale * Eigen::Vector3d(1e-5, -8e-6, 6e-6);
      Vector sigma(8);
      sigma << 2, 2, 2, 1e-4, 1e-4, 1e-4, ideal ? 0 : 1e-5,
          ideal ? 0 : options.clockFrequencySigma;
      filter.setInitialState(initial,
                             sigma.array().square().matrix().asDiagonal(), t0);
      auto prior = csv(folder / "initial_prior.csv");
      prior << "parameter,mean,standard_deviation,truth\n";
      const std::array<std::string, 8> parameterNames{
          "rx_km",   "ry_km",   "rz_km",        "vx_km_s",
          "vy_km_s", "vz_km_s", "clock_bias_s", "fractional_frequency"};
      Vector trueInitial(8);
      trueInitial.head<6>() = truth.relative(t0);
      trueInitial[6] = clockTruth.nodes.front().bias_s;
      trueInitial[7] = clockTruth.nodes.front().fractional_frequency;
      for (int i = 0; i < 8; ++i)
        prior << parameterNames[i] << ',' << initial[i] << ',' << sigma[i]
              << ',' << trueInitial[i] << '\n';
      auto clockFile = csv(folder / "clock_truth.csv");
      clockFile << "elapsed_s,bias_s,fractional_frequency\n";
      for (std::size_t i = 0; i < clockTruth.nodes.size(); ++i)
        clockFile << i * 60 << ',' << clockTruth.nodes[i].bias_s << ','
                  << clockTruth.nodes[i].fractional_frequency << '\n';
      auto trajectory = csv(folder / "trajectory.csv"),
           innovations = csv(folder / "innovations.csv"),
           covariances = csv(folder / "covariance.csv");
      trajectory
          << "elapsed_s,truth_tag_elapsed_s,tracking,x_km,y_km,z_km,vx_km_s,vy_"
             "km_s,vz_km_s,truth_x_km,truth_y_km,truth_z_km,truth_vx_km_s,"
             "truth_vy_km_s,truth_vz_km_s,radial_error_m,along_error_m,cross_"
             "error_m,radial_velocity_error_mm_s,along_velocity_error_mm_s,"
             "cross_velocity_error_mm_s,position_error_m,velocity_error_mm_s,"
             "formal_position_sigma_m,formal_velocity_sigma_mm_s,clock_bias_"
             "error_s,clock_frequency_error,clock_bias_sigma_s,clock_frequency_"
             "sigma,clock_bias_s,clock_frequency,truth_clock_bias_s,truth_"
             "clock_frequency,orbit_nees,joint_nees\n";
      innovations
          << "elapsed_s,station,range_active,range_innovation_m,doppler_"
             "innovation_hz,range_innovation_sigma_m,doppler_innovation_sigma_"
             "hz,range_whitened,doppler_whitened,nis,dimension,accepted\n";
      covariances << "elapsed_s";
      for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
          covariances << ",P" << i << j;
      covariances << '\n';
      double positionSum = 0, velocitySum = 0, nisSum = 0, lastPosition = 0,
             lastVelocity = 0;
      double maxStencilFraction = 0;
      std::size_t samples = 0, observationIndex = 0, updates = 0;
      const auto exportState = [&](double tagElapsed, bool tracking) {
        const double t = filter.epoch();
        const auto x = filter.state();
        const auto &p = filter.covariance();
        const auto exact = truth.relative(t);
        const auto c = clockTruth.stateAt(t);
        Vector error = x;
        error.head<6>() -= exact;
        error[6] -= c.bias_s;
        error[7] -= c.fractional_frequency;
        const Eigen::Vector3d
            radial = exact.head<3>().normalized(),
            cross = exact.head<3>().cross(exact.tail<3>()).normalized(),
            along = cross.cross(radial);
        Eigen::Matrix3d rotation;
        rotation.row(0) = radial.transpose();
        rotation.row(1) = along.transpose();
        rotation.row(2) = cross.transpose();
        const Eigen::Vector3d dr = rotation * error.head<3>() * 1000,
                              dv = rotation * error.segment<3>(3) * 1e6;
        lastPosition = dr.norm();
        lastVelocity = dv.norm();
        positionSum += dr.squaredNorm();
        velocitySum += dv.squaredNorm();
        ++samples;
        Vector scales(8);
        scales << 2, 2, 2, 1e-4, 1e-4, 1e-4, 1e-5, 1e-10;
        const auto inv = scales.cwiseInverse().eval();
        const Eigen::MatrixXd normalized =
            inv.asDiagonal() * p * inv.asDiagonal();
        const Vector normalizedError = inv.asDiagonal() * error;
        const double orbitNees = normalizedError.head<6>().dot(
            normalized.topLeftCorner<6, 6>().llt().solve(
                normalizedError.head<6>()));
        const double jointNees =
            ideal
                ? orbitNees
                : normalizedError.dot(normalized.llt().solve(normalizedError));
        trajectory << t - t0 << ',' << tagElapsed << ',' << tracking;
        for (int i = 0; i < 6; ++i)
          trajectory << ',' << x[i];
        for (int i = 0; i < 6; ++i)
          trajectory << ',' << exact[i];
        for (int i = 0; i < 3; ++i)
          trajectory << ',' << dr[i];
        for (int i = 0; i < 3; ++i)
          trajectory << ',' << dv[i];
        trajectory << ',' << lastPosition << ',' << lastVelocity << ','
                   << std::sqrt(p.topLeftCorner<3, 3>().trace()) * 1000 << ','
                   << std::sqrt(p.block<3, 3>(3, 3).trace()) * 1e6 << ','
                   << error[6] << ',' << error[7] << ',' << std::sqrt(p(6, 6))
                   << ',' << std::sqrt(p(7, 7)) << ',' << x[6] << ',' << x[7]
                   << ',' << c.bias_s << ',' << c.fractional_frequency << ','
                   << orbitNees << ',' << jointNees << '\n';
        covariances << t - t0;
        for (int i = 0; i < 8; ++i)
          for (int j = 0; j < 8; ++j)
            covariances << ',' << p(i, j);
        covariances << '\n';
      };
      exportState(0, false);
      std::cout << kind << ": " << observations.size() << " Doppler, "
                << rangeCount << " range; running causal EKF." << std::endl;
      for (int tick = 1; tick < static_cast<int>(duration / 60); ++tick) {
        const double elapsed = tick * 60.;
        const bool tracking = observationIndex < observations.size() &&
                              observations[observationIndex].elapsed == elapsed;
        if (!tracking && tick % 10 != 0 &&
            tick != static_cast<int>(duration / 60) - 1)
          continue;
        const long double tag = tracking ? observations[observationIndex].tag
                                         : receiverTruth.reading(t0 + elapsed);
        const double currentEpoch = filter.epoch();
        const Vector previous = filter.state();
        const radio::EndpointClock estimatedClock{
            t0,
            [&, previous, currentEpoch](double t) {
              const double dt = t - currentEpoch;
              return previous[6] + previous[7] * dt +
                     .5 * parameters.frequency_drift_per_s * dt * dt;
            },
            nominalOffset};
        // Epoch assignment uses the raw clock tag and current estimate only.
        const double epoch = estimatedClock.coordinateEpoch(tag);
        filter.predictTo(epoch, propagate);
        if (tracking) {
          const auto &a = observations[observationIndex++];
          const Vector before = filter.state();
          const cassini_benchmark::ReceptionTrajectory local(
              epoch, before, acceleration, saturn, options.maxStep);
          Vector steps(8); // Large symmetric stencils on affine trajectories
                           // resolve GHz count arithmetic.
          steps << 1000, 1000, 1000, 1e-2, 1e-2, 1e-2, 1e-2, 1e-5;
          const radio::OneWayMeasurement measurement(
              links[a.station], *stations[a.station], tag, 60, t0, ramp,
              ramp.frequencyHz, ground, nominalOffset,
              parameters.frequency_drift_per_s, local.factory(), steps);
          const auto prediction = measurement(epoch, before);
          if (updates % 500 == 0) {
            const radio::OneWayMeasurement half(
                links[a.station], *stations[a.station], tag, 60, t0, ramp,
                ramp.frequencyHz, ground, nominalOffset,
                parameters.frequency_drift_per_s, local.factory(), steps / 2);
            const auto alternate = half(epoch, before);
            const auto disagreement =
                (prediction.jacobian - alternate.jacobian) *
                filter.covariance() *
                (prediction.jacobian - alternate.jacobian).transpose();
            maxStencilFraction =
                std::max(maxStencilFraction,
                         std::max(std::sqrt(std::max(0., disagreement(0, 0))) /
                                      (options.rangeSigmaM / 1000),
                                  std::sqrt(std::max(0., disagreement(1, 1))) /
                                      options.dopplerSigmaHz));
          }
          const int m = a.range ? 2 : 1;
          Eigen::MatrixXd r = Eigen::MatrixXd::Zero(m, m);
          Vector observed(m);
          if (a.range) {
            r(0, 0) = std::pow(options.rangeSigmaM / 1000, 2);
            r(1, 1) = std::pow(options.dopplerSigmaHz, 2);
            observed = a.value;
          } else {
            r(0, 0) = std::pow(options.dopplerSigmaHz, 2);
            observed[0] = a.value[1];
          }
          const auto diagnostic =
              filter.update(epoch, observed, r, [&](double, const Vector &) {
                if (a.range)
                  return prediction;
                return filters::MeasurementPrediction{
                    prediction.value.tail(1),
                    prediction.jacobian.bottomRows(1)};
              });
          innovations << epoch - t0 << ',' << names[a.station] << ',' << a.range
                      << ',';
          if (a.range)
            innovations << diagnostic.innovation[0] * 1000;
          innovations << ',' << diagnostic.innovation[m - 1] << ',';
          if (a.range)
            innovations << std::sqrt(diagnostic.covariance(0, 0)) * 1000;
          innovations << ',' << std::sqrt(diagnostic.covariance(m - 1, m - 1))
                      << ',';
          if (a.range)
            innovations << diagnostic.whitenedInnovation[0];
          innovations << ',' << diagnostic.whitenedInnovation[m - 1] << ','
                      << diagnostic.normalizedInnovationSquared << ',' << m
                      << ',' << diagnostic.accepted << '\n';
          nisSum += diagnostic.normalizedInnovationSquared;
          ++updates;
          if (updates % 1000 == 0)
            std::cout << kind << ": " << updates << " updates, position error "
                      << (filter.state().head<3>() -
                          truth.relative(epoch).head<3>())
                                 .norm() *
                             1000
                      << " m" << std::endl;
        }
        // Fixed 600 s diagnostic grid, including tracking gaps (not a truth
        // input).
        if (tick % 10 == 0 || tick == static_cast<int>(duration / 60) - 1)
          exportState(elapsed, tracking);
      }
      if (updates != observations.size() || maxStencilFraction > .01)
        throw std::runtime_error(
            "Observation/Jacobian validation failed: updates=" +
            std::to_string(updates) +
            " observations=" + std::to_string(observations.size()) +
            " stencil_fraction=" + std::to_string(maxStencilFraction));
      auto summary = csv(folder / "summary.json");
      summary << "{\n\"clock_case\":\"" << kind
              << "\",\n\"duration_days\":" << options.days
              << ",\n\"seed\":" << options.seed
              << ",\n\"noise_enabled\":" << (options.noise ? "true" : "false")
              << ",\n\"estimated_parameters\":8,\n\"truth_source\":\"results/"
                 "cassini_two_way_wls/pass_schedule/truth/"
                 "trajectory.csv:multibody_j2_srp_thermal\",\n\"clock_"
                 "configuration\":\""
              << clock.configurationName() << "\",\n\"clock_drift_per_s\":"
              << parameters.frequency_drift_per_s
              << ",\n\"clock_q_bias_s\":" << parameters.q_bias_s
              << ",\n\"clock_q_frequency_per_s\":"
              << parameters.q_frequency_per_s
              << ",\n\"range_sigma_m\":" << options.rangeSigmaM
              << ",\n\"doppler_sigma_hz\":" << options.dopplerSigmaHz
              << ",\n\"initial_guess_scale\":" << options.guessScale
              << ",\n\"clock_frequency_prior_sigma\":"
              << (ideal ? 0 : options.clockFrequencySigma)
              << ",\n\"max_step_s\":" << options.maxStep
              << ",\n\"doppler_observations\":" << updates
              << ",\n\"range_observations\":" << rangeCount
              << ",\n\"visibility_rejected_candidates\":" << rejected
              << ",\n\"diagnostic_samples\":" << samples
              << ",\n\"position_rms_m\":" << std::sqrt(positionSum / samples)
              << ",\n\"velocity_rms_mm_s\":" << std::sqrt(velocitySum / samples)
              << ",\n\"endpoint_position_error_m\":" << lastPosition
              << ",\n\"endpoint_velocity_error_mm_s\":" << lastVelocity
              << ",\n\"mean_nis\":" << nisSum / updates
              << ",\n\"maximum_half_stencil_noise_fraction\":"
              << maxStencilFraction << ",\n\"wall_seconds\":"
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - wall)
                     .count()
              << "\n}\n";
      std::cout << kind << ": RMS " << std::sqrt(positionSum / samples)
                << " m, endpoint " << lastPosition << " m, stencil fraction "
                << maxStencilFraction << std::endl;
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
