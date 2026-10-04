#include "Clocks/ClockTruthSimulator.hpp"
#include "Clocks/DSAC.hpp"
#include "dynamics/SpiceEarthEnvironment.hpp"
#include "filters/CartesianPropagator.hpp"
#include "filters/EKF.hpp"
#include "observations/GeometricRadiometricModel.hpp"
#include "perturbations/Eclipse.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/SRP.hpp"
#include "stations/StationCatalog.hpp"
#include <Eigen/Eigenvalues>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>

using namespace fd::filters;
namespace {
constexpr double step = 60, duration = 8 * 3600;

AccelerationFunction makeForces(const fd::dynamics::SpiceEarthEnvironment& environment,
                                double startEpoch, bool twoBodyOnly) {
    auto gravity = fd::perturbations::pointMassGravity(environment.earthMu());
    if (twoBodyOnly)
        return gravity;
    const fd::perturbations::SolarRadiationPressure srp("EARTH", "J2000", 1.3, 20, 1000);
    return [&, startEpoch, gravity, srp](double elapsed, const Eigen::Vector3d& position,
                                         const Eigen::Vector3d& velocity) {
        const auto sun = environment.sunPosition(startEpoch + elapsed);
        auto force = gravity(elapsed, position, velocity);
        const auto tidal = fd::perturbations::thirdBodyGravity(environment.sunMu(), position, sun);
        const auto radiation = srp.evaluateWithShadow(position, sun, environment.earthRadius(),
                                                      environment.sunRadius());
        force.acceleration += tidal.acceleration + radiation.acceleration;
        force.positionJacobian += tidal.positionJacobian + radiation.positionJacobian;
        return force;
    };
}

void run(int n, const std::filesystem::path& output, bool gaps, bool twoBodyOnly,
         const fd::dynamics::SpiceEarthEnvironment& environment, double startEpoch,
         const std::vector<od::Station>& stations) {
    constexpr double radius = 20000, rad = std::numbers::pi / 180;
    const auto clock = fd::clocks::DSAC::shortTermWhiteFmBaseline();
    fd::clocks::ClockTruthSimulator clockTruth(clock, 2027);
    CartesianPropagationConfig config;
    config.integrator.absoluteTolerance = 1e-11;
    config.integrator.relativeTolerance = 1e-11;
    config.integrator.initialStep = 30;
    config.integrator.maximumStep = 10; // Resolve the short penumbra transitions.
    config.clock = clock.parameters();
    const CartesianPropagator propagate(makeForces(environment, startEpoch, twoBodyOnly), config);
    EKF filter(n == 6 ? StateLayout::Orbit : StateLayout::OrbitClock);
    Eigen::VectorXd truth = Eigen::VectorXd::Zero(n);
    const double speed = std::sqrt(environment.earthMu() / radius);
    truth.head<6>() << radius, 0, 0, 0, speed * std::cos(25 * rad), speed * std::sin(25 * rad);
    if (n == 8) {
        truth[6] = 2e-7;
        truth[7] = 2e-11;
    }
    Eigen::VectorXd initial = truth;
    initial.head<6>() +=
        (Eigen::Matrix<double, 6, 1>() << .6, -.4, .3, 1e-4, -8e-5, 6e-5).finished();
    Eigen::VectorXd sigma = Eigen::VectorXd::Constant(n, 1e-3);
    sigma.head<3>().setConstant(1);
    if (n == 8) {
        initial.tail<2>().setZero();
        sigma[6] = 1e-6;
        sigma[7] = 1e-10;
    }
    filter.setInitialState(initial, sigma.array().square().matrix().asDiagonal(), 0);
    const Eigen::Vector2d measurementSigma(.01, 1e-7);
    const Eigen::Matrix2d measurementCovariance =
        measurementSigma.array().square().matrix().asDiagonal();
    std::mt19937_64 measurementGenerator(2026);
    std::normal_distribution<double> measurementNormal;
    std::ofstream csv(output / ((std::string("ekf_") + std::to_string(n)) +
                                (gaps ? "_gaps.csv" : "_continuous.csv")));
    if (!csv)
        throw std::runtime_error("Cannot create EKF diagnostics CSV.");
    csv << std::setprecision(17)
        << "time_s,tracking,station,illumination_truth,eclipse_truth,"
           "illumination_estimate,eclipse_estimate,position_error_km,position_"
           "bound_3sigma_km";
    for (int j = 0; j < n; ++j)
        csv << ",error_" << j << ",sigma_" << j;
    csv << ",range_innovation_km,rate_innovation_km_s,range_innovation_sigma_km,"
           "rate_innovation_"
           "sigma_km_s,white_range,white_rate,nis\n";
    std::size_t updates = 0;
    for (int k = 0; k <= static_cast<int>(duration / step); ++k) {
        const double time = k * step;
        if (k > 0) {
            // Propagate orbital truth separately so the clock simulator owns every
            // clock increment.
            truth.head<6>() = propagate(time - step, time, truth.head<6>()).state;
            if (n == 8) {
                const auto next = clockTruth.step({truth[6], truth[7]}, step);
                truth[6] = next.bias_s;
                truth[7] = next.fractional_frequency;
            }
            filter.predictTo(time, propagate);
        }
        const bool scheduled = !gaps || (!(time >= 2 * 3600 && time < 4 * 3600) &&
                                         !(time >= 6 * 3600 && time < 7 * 3600));
        // Draw on every epoch so contact schedules share the same noise
        // realizations.
        const Eigen::Vector2d measurementNoise(
            measurementSigma[0] * measurementNormal(measurementGenerator),
            measurementSigma[1] * measurementNormal(measurementGenerator));
        int selected = -1;
        double bestElevation = std::sin(10 * rad);
        fd::dynamics::CartesianState station;
        for (std::size_t i = 0; i < stations.size(); ++i) {
            const auto candidate = environment.stationState(stations[i], startEpoch + time);
            const Eigen::Vector3d line = truth.head<3>() - candidate.positionKm;
            const double elevation =
                line.normalized().dot(environment.stationUp(stations[i], startEpoch + time));
            if (elevation > bestElevation) {
                bestElevation = elevation;
                selected = i;
                station = candidate;
            }
        }
        const bool tracking = scheduled && selected >= 0;
        InnovationDiagnostics diagnostics;
        if (tracking) {
            const auto model = [station](double, const Eigen::VectorXd& x) {
                return fd::observations::geometricRadiometricPrediction(x, station);
            };
            Eigen::Vector2d observed = model(time, truth).value;
            observed += measurementNoise;
            diagnostics = filter.update(time, observed, measurementCovariance, model);
            ++updates;
        }
        const Eigen::VectorXd error = filter.state() - truth;
        const Eigen::Matrix3d positionCovariance = filter.covariance().topLeftCorner<3, 3>();
        const double bound =
            3 * std::sqrt(Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(positionCovariance)
                              .eigenvalues()
                              .maxCoeff());
        const auto sun = environment.sunPosition(startEpoch + time);
        const auto trueShadow = fd::perturbations::evaluateEclipse(
            truth.head<3>(), sun, environment.earthRadius(), environment.sunRadius());
        const auto estimatedShadow = fd::perturbations::evaluateEclipse(
            filter.state().head<3>(), sun, environment.earthRadius(), environment.sunRadius());
        csv << time << ',' << tracking << ',' << (tracking ? selected : -1) << ','
            << trueShadow.illumination << ',' << static_cast<int>(trueShadow.flag) << ','
            << estimatedShadow.illumination << ',' << static_cast<int>(estimatedShadow.flag) << ','
            << error.head<3>().norm() << ',' << bound;
        for (int j = 0; j < n; ++j)
            csv << ',' << error[j] << ',' << std::sqrt(filter.covariance()(j, j));
        if (tracking)
            csv << ',' << diagnostics.innovation[0] << ',' << diagnostics.innovation[1] << ','
                << std::sqrt(diagnostics.covariance(0, 0)) << ','
                << std::sqrt(diagnostics.covariance(1, 1)) << ','
                << diagnostics.whitenedInnovation[0] << ',' << diagnostics.whitenedInnovation[1]
                << ',' << diagnostics.normalizedInnovationSquared;
        else
            csv << ",nan,nan,nan,nan,nan,nan,nan";
        csv << '\n';
    }
    if (!csv)
        throw std::runtime_error("Cannot write EKF diagnostics CSV.");
    if (updates == 0)
        throw std::runtime_error("EKF demo has no visible tracking observations.");
    std::cout << n << " states, " << (gaps ? "gaps" : "continuous") << ": " << updates
              << " updates; final position error "
              << (filter.state() - truth).head<3>().norm() * 1000 << " m\n";
}
} // namespace

int main(int argc, char** argv) {
    try {
        std::filesystem::path output = "Output EKF";
        int selected = 0;
        bool twoBodyOnly = false;
        std::filesystem::path kernelDirectory =
            std::filesystem::path(DEEPNAV_SOURCE_DIR) / "Kernels";
        std::string startUtc = "2024-03-20T00:00:00";
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument == "--two-body")
                twoBodyOnly = true;
            else if (argument == "--kernels" && i + 1 < argc)
                kernelDirectory = argv[++i];
            else if (argument == "--start-utc" && i + 1 < argc)
                startUtc = argv[++i];
            else if (argument == "--output" && i + 1 < argc)
                output = argv[++i];
            else if (argument == "--states" && i + 1 < argc) {
                const std::string value = argv[++i];
                if (value != "6" && value != "8")
                    throw std::invalid_argument("--states must be 6 or 8.");
                selected = value == "6" ? 6 : 8;
            } else
                throw std::invalid_argument(
                    "Usage: ekf_demo [--states 6|8] [--output DIR] [--two-body] "
                    "[--kernels DIR] [--start-utc UTC]");
        }
        const fd::dynamics::SpiceEarthEnvironment environment(kernelDirectory);
        const double startEpoch = environment.epochTdb(startUtc);
        std::vector<od::Station> stations;
        for (const char* name : {"DSS-43", "DSS-63", "DSS-14"}) {
            stations.push_back(
                od::buildStationFromKernel(name, od::stationNaifIdFromName(name), startEpoch));
            (void)environment.stationState(stations.back(), startEpoch + duration);
            (void)environment.stationUp(stations.back(), startEpoch + duration);
        }
        (void)environment.sunPosition(startEpoch + duration);
        std::filesystem::create_directories(output);
        std::ofstream metadata(output / "scenario.json");
        metadata << std::setprecision(17)
                 << "{\n  \"integrator\": \"RKF45\",\n  \"force_model\": \""
                 << (twoBodyOnly ? "Earth point mass"
                                 : "Earth + Sun gravity + eclipse-modulated SRP")
                 << "\",\n  \"start_epoch_tdb\": " << startEpoch
                 << ",\n  \"sun_geometry\": \"CSPICE DE442, geometric J2000\",\n"
                 << "  \"stations\": [\"DSS-43\", \"DSS-63\", \"DSS-14\"],\n"
                 << "  \"eclipse_model\": \"Spherical angular-disk overlap, Earth "
                    "equatorial radius\",\n"
                 << "  \"clock_model\": \"DSAC_inspired_short_term_white_FM\",\n"
                 << "  \"clock_truth_seed\": 2027,\n  \"measurement_seed\": 2026,\n"
                 << "  \"srp_cr\": 1.3,\n  \"srp_area_m2\": 20,\n  "
                    "\"srp_mass_kg\": 1000\n}\n";
        if (!metadata)
            throw std::runtime_error("Cannot write EKF scenario metadata.");
        std::cout << "RKF45; " << (twoBodyOnly ? "Earth point mass" : "Earth + solar gravity + SRP")
                  << '\n';
        for (int n : {6, 8})
            if (selected == 0 || selected == n) {
                run(n, output, false, twoBodyOnly, environment, startEpoch, stations);
                run(n, output, true, twoBodyOnly, environment, startEpoch, stations);
            }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
