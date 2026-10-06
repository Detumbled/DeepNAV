#include "dynamics/IntegratedStateProvider.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "filters/CartesianPropagator.hpp"
#include "filters/EKF.hpp"
#include "observations/radiometric/PropagationCorrections.hpp"
#include "observations/radiometric/SpiceLightTime.hpp"
#include "observations/radiometric/TwoWayMeasurement.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/J2.hpp"
#include "simulation/CassiniRadioKernels.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <random>

namespace radio = fd::observations::radiometric;
namespace {
using State6 = Eigen::Matrix<double, 6, 1>;
struct Options {
    double hours{24}, step{600}, count{60}, noiseDensity{1e-5};
    std::string truthMode{"spice"};
    std::string utc{"2004-10-10T18:13:06.046"};
    std::filesystem::path output{"Output EKF/Cassini two-way"};
};
Options parse(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--help") {
            std::cout << "cassini_two_way_ekf_demo [--hours H] [--step-s S] [--count-s S] "
                         "[--start-utc UTC] [--output DIR] [--truth spice|matched]\n"
                         "  [--acceleration-noise-density M/S^(3/2)]\n";
            std::exit(0);
        }
        if (++i >= argc)
            throw std::invalid_argument("Missing value for " + key);
        if (key == "--start-utc")
            o.utc = argv[i];
        else if (key == "--output")
            o.output = argv[i];
        else if (key == "--truth")
            o.truthMode = argv[i];
        else {
            std::size_t used;
            const std::string text = argv[i];
            const double x = std::stod(text, &used);
            if (used != text.size() || !std::isfinite(x))
                throw std::invalid_argument("Invalid value");
            if (key == "--hours")
                o.hours = x;
            else if (key == "--step-s")
                o.step = x;
            else if (key == "--count-s")
                o.count = x;
            else if (key == "--acceleration-noise-density")
                o.noiseDensity = x;
            else
                throw std::invalid_argument("Unknown option " + key);
        }
    }
    if (o.hours <= 0 || o.hours > 72 || o.step <= 0 || o.count <= 0 || o.count > o.step ||
        o.count > 300 || o.noiseDensity < 0 || (o.truthMode != "spice" && o.truthMode != "matched"))
        throw std::invalid_argument("Positive arc <=72h; count <= step and <=300s required.");
    return o;
}
double gm(const char *body) {
    SpiceInt n;
    double value;
    bodvrd_c(body, "GM", 1, &n, &value);
    od::throwIfSpiceFailed(std::string("GM ") + body);
    return value;
}
State6 relativeSpice(const char *body, double epoch) {
    double x[6], lt;
    spkezr_c(body, epoch, "J2000", "NONE", "SATURN", x, &lt);
    od::throwIfSpiceFailed(std::string("Saturn-relative state ") + body);
    return Eigen::Map<State6>(x);
}
std::pair<double, double> elevations(const std::string &station, const radio::TwoWaySolution &s) {
    return {radio::spiceStationElevationRadians(station, s.uplink.emissionTdb,
                                                s.uplink.receiver.positionKm -
                                                    s.uplink.emitter.positionKm),
            radio::spiceStationElevationRadians(station, s.downlink.receptionTdb,
                                                s.downlink.emitter.positionKm -
                                                    s.downlink.receiver.positionKm)};
}
} // namespace
int main(int argc, char **argv) {
    try {
        const auto o = parse(argc, argv);
        const auto wallStart = std::chrono::steady_clock::now();
        od::SpiceErrorModeGuard errors("RETURN", "NONE");
        fd::dynamics::SpiceKernelSet kernels(fd::simulation::cassiniRadioKernels(
            std::filesystem::path(DEEPNAV_SOURCE_DIR) / "Kernels"));
        const radio::SpiceTdbClock clock;
        double start;
        str2et_c(o.utc.c_str(), &start);
        od::throwIfSpiceFailed("Arc start UTC");
        const double duration = o.hours * 3600, lookback = 5000, initialEpoch = start - lookback;
        fd::dynamics::SpiceStateProvider saturn("SATURN"), sun("SUN"), cassini("CASSINI");
        const double saturnGm = gm("SATURN"), solarGm = gm("SUN");
        std::vector<std::pair<const char *, double>> perturbingBodies;
        for (const char *body : {"SUN", "MIMAS", "ENCELADUS", "TETHYS", "DIONE", "RHEA", "TITAN",
                                 "HYPERION", "IAPETUS"})
            perturbingBodies.emplace_back(body, gm(body));
        double radii[3];
        SpiceInt n;
        bodvrd_c("SATURN", "RADII", 3, &n, radii);
        double rotation[3][3];
        pxform_c("IAU_SATURN", "J2000", initialEpoch, rotation);
        od::throwIfSpiceFailed("Saturn radius and pole");
        const Eigen::Vector3d pole(rotation[0][2], rotation[1][2], rotation[2][2]);
        const auto central = fd::perturbations::pointMassGravity(saturnGm);
        // Reduced estimator dynamics; not a reconstruction of flight OD forces.
        // J2 is the rounded NASA fact-sheet value.
        const fd::perturbations::AccelerationFunction force =
            [&](double t, const Eigen::Vector3d &r, const Eigen::Vector3d &v) {
                auto total = central(t, r, v);
                for (const auto &body : perturbingBodies) {
                    const auto third = fd::perturbations::thirdBodyGravity(
                        body.second, r, relativeSpice(body.first, t).head<3>());
                    total.acceleration += third.acceleration;
                    total.positionJacobian += third.positionJacobian;
                }
                const auto j2 = fd::perturbations::j2Gravity(saturnGm, radii[0], .016298, r, pole);
                total.acceleration += j2.acceleration;
                total.positionJacobian += j2.positionJacobian;
                return total;
            };
        od::RKF45Integrator::Options integration;
        integration.absoluteTolerance = 1e-10;
        integration.relativeTolerance = 1e-12;
        integration.initialStep = 30;
        integration.maximumStep = 60;
        const State6 truthInitial = relativeSpice("CASSINI", initialEpoch);
        std::unique_ptr<fd::dynamics::IntegratedStateProvider> matched;
        if (o.truthMode == "matched")
            matched = std::make_unique<fd::dynamics::IntegratedStateProvider>(
                initialEpoch, truthInitial, start + duration + o.step, force, integration, "J2000",
                "SSB", &saturn);
        const fd::dynamics::StateProvider &truth =
            matched ? static_cast<const fd::dynamics::StateProvider &>(*matched) : cassini;
        fd::filters::CartesianPropagationConfig config;
        config.integrator = integration;
        // SNC models residual acceleration uncertainty separately from R.
        config.accelerationDiffusion =
            std::pow(o.noiseDensity / 1000, 2) * Eigen::Matrix3d::Identity();
        const fd::filters::CartesianPropagator propagate(
            [&](double elapsed, const Eigen::Vector3d &r, const Eigen::Vector3d &v) {
                return force(initialEpoch + elapsed, r, v);
            },
            config);
        fd::filters::EKF filter;
        State6 initial = truthInitial;
        initial += (State6() << 1, -.8, .5, 1e-5, -8e-6, 6e-6).finished();
        State6 sigma;
        sigma << 2, 2, 2, 1e-4, 1e-4, 1e-4;
        filter.setInitialState(initial, sigma.array().square().matrix().asDiagonal(), 0);
        const radio::LinkConfig transponder{880.0 / 749.0, 1e-6};
        const radio::TwoWayLink geometric(radio::libraryReceptionSolver(), transponder);
        const radio::FrequencyRamp ramp{start, 7.175e9, 0};
        const double referenceHz = transponder.turnaroundRatio * ramp.frequencyHz;
        const Eigen::Vector2d noiseSigma(.01, .01); // 10 m range; 0.01 Hz counted Doppler.
        const Eigen::Matrix2d measurementCovariance =
            noiseSigma.array().square().matrix().asDiagonal();
        State6 differentiation;
        differentiation << .1, .1, .1, 1e-5, 1e-5, 1e-5;
        std::mt19937_64 generator(2026);
        std::normal_distribution<double> normal;
        std::filesystem::create_directories(o.output);
        std::ofstream csv(o.output / "diagnostics.csv");
        csv.exceptions(std::ios::badbit | std::ios::failbit);
        csv << std::setprecision(17);
        csv << "reception_elapsed_s,state_elapsed_s,station,tracking,position_error_m,"
               "running_position_rms_m,nis,range_innovation_m,doppler_innovation_hz\n";
        double sumSquared = 0, nisSum = 0, firstError = 0, lastError = 0, maxNis = 0;
        std::size_t updates = 0, epochs = 0;
        std::cout
            << "Cassini two-way EKF (truth=" << o.truthMode << "): " << o.hours << " h, cadence "
            << o.step << " s, count " << o.count
            << " s. Saturn + J2 + solar/eight-moon tides; SNC density " << o.noiseDensity
            << " m/s^(3/2).\n"
               "6-state Saturn-relative EKF; SSB link geometry; fixed 5000 s state lookback.\n";
        const int samples = static_cast<int>(std::floor(duration / o.step)) + 1;
        for (int k = 0; k < samples; ++k) {
            const double elapsed = k * o.step, reception = start + elapsed,
                         stateEpoch = reception - lookback;
            filter.predictTo(elapsed, propagate);
            std::string chosen;
            double best = -1e9;
            for (const std::string name : {"DSS-14", "DSS-43", "DSS-63"}) {
                fd::dynamics::SpiceStateProvider station(name);
                const auto a =
                    elevations(name, geometric.solve(reception - o.count, station, truth));
                const auto b = elevations(name, geometric.solve(reception, station, truth));
                const double minimum = std::min({a.first, a.second, b.first, b.second});
                if (minimum > best) {
                    best = minimum;
                    chosen = name;
                }
            }
            double nis = std::numeric_limits<double>::quiet_NaN();
            Eigen::Vector2d innovation = Eigen::Vector2d::Constant(nis);
            const bool tracking = best >= 10 * std::numbers::pi / 180;
            if (tracking) {
                fd::dynamics::SpiceStateProvider station(chosen);
                const auto correction = [&](const radio::LegSolution &leg, const auto &emitter,
                                            const auto &) {
                    const bool up = &emitter == &station;
                    const auto direction = up ? leg.receiver.positionKm - leg.emitter.positionKm
                                              : leg.emitter.positionKm - leg.receiver.positionKm;
                    const double elevation = radio::spiceStationElevationRadians(
                        chosen, up ? leg.emissionTdb : leg.receptionTdb, direction);
                    if (elevation <= 0)
                        throw std::runtime_error("Corrected link below horizon");
                    return radio::pointMassShapiroSeconds(leg, sun, solarGm) +
                           radio::zenithPathDelaySeconds(2.3, elevation);
                };
                const radio::TwoWayLink link(radio::libraryReceptionSolver({}, correction),
                                             transponder);
                const auto truthCount =
                    link.count(reception, o.count, station, truth, ramp, referenceHz, clock);
                Eigen::Vector2d observed(truthCount.end.rangeKm(), truthCount.residualHz);
                observed[0] += noiseSigma[0] * normal(generator);
                observed[1] += noiseSigma[1] * normal(generator);
                const radio::TrajectoryFactory factory = [&](double epoch,
                                                             const Eigen::VectorXd &x) {
                    return std::make_unique<fd::dynamics::IntegratedStateProvider>(
                        epoch, State6(x), reception + 1, force, integration, "J2000", "SSB",
                        &saturn);
                };
                const radio::TwoWayMeasurement measurement(link, station, reception, o.count, ramp,
                                                           referenceHz, clock, factory,
                                                           differentiation);
                // Filter epoch is a driver coordinate, not Earth reception. The
                // trajectory factory maps it to TDB and evaluates every retarded leg.
                const auto diagnostics = filter.update(elapsed, observed, measurementCovariance,
                                                       [&](double t, const Eigen::VectorXd &x) {
                                                           return measurement(initialEpoch + t, x);
                                                       });
                nis = diagnostics.normalizedInnovationSquared;
                innovation = diagnostics.innovation;
                nisSum += nis;
                maxNis = std::max(maxNis, nis);
                ++updates;
            }
            const Eigen::Vector3d actual = matched ? matched->relativeStateAt(stateEpoch).positionKm
                                                   : relativeSpice("CASSINI", stateEpoch).head<3>();
            lastError = (filter.state().head<3>() - actual).norm() * 1000;
            if (k == 0)
                firstError = lastError;
            sumSquared += lastError * lastError;
            ++epochs;
            const double rms = std::sqrt(sumSquared / epochs);
            csv << elapsed << ',' << stateEpoch - initialEpoch << ',' << chosen << ',' << tracking
                << ',' << lastError << ',' << rms << ',' << nis << ',' << innovation[0] * 1000
                << ',' << innovation[1] << '\n';
            if (k % 12 == 0 || k == samples - 1)
                std::cout << std::fixed << std::setprecision(3) << "t=" << elapsed / 3600
                          << " h, updates=" << updates << ", position error=" << lastError
                          << " m, running 3D position RMS=" << rms << " m, NIS=" << nis << '\n';
        }
        csv.close();
        if (!updates)
            throw std::runtime_error("No visible two-way updates on this arc.");
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - wallStart).count();
        std::ofstream summary(o.output / "summary.json");
        summary << std::setprecision(17) << "{\n  \"arc_hours\": " << o.hours
                << ",\n  \"epochs\": " << epochs << ",\n  \"updates\": " << updates
                << ",\n  \"measurement_dimension\": 2"
                << ",\n  \"initial_postfit_position_error_m\": " << firstError
                << ",\n  \"final_position_error_m\": " << lastError
                << ",\n  \"running_position_rms_m\": " << std::sqrt(sumSquared / epochs)
                << ",\n  \"mean_nis\": " << nisSum / updates << ",\n  \"max_nis\": " << maxNis
                << ",\n  \"runtime_s\": " << seconds << ",\n  \"truth\": \"" << o.truthMode << "\""
                << ",\n  \"acceleration_noise_density_m_s32\": " << o.noiseDensity
                << ",\n  \"noise_seed\": 2026"
                << ",\n  \"range_sigma_m\": 10,\n  \"doppler_sigma_hz\": 0.01\n}\n";
        std::cout << "Finished " << epochs << " epochs, " << updates << " updates (" << 2 * updates
                  << " scalar observations). Mean NIS=" << nisSum / updates
                  << "; runtime=" << seconds
                  << " s. Diagnostics: " << (o.output / "diagnostics.csv") << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Cassini two-way EKF: " << e.what() << '\n';
        return 1;
    }
}
