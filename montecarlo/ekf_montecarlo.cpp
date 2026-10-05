#include "Clocks/ClockTruthSimulator.hpp"
#include "Clocks/DSAC.hpp"
#include "Statistics.hpp"
#include "dynamics/SpiceEarthEnvironment.hpp"
#include "filters/CartesianPropagator.hpp"
#include "filters/EKF.hpp"
#include "observations/GeometricRadiometricModel.hpp"
#include "simulation/EarthOrbitScenario.hpp"
#include "stations/StationCatalog.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>

using namespace fd::filters;
namespace {
constexpr double step = 60, duration = 8 * 3600, srpSigma = .2;
constexpr int epochs = static_cast<int>(duration / step) + 1;
constexpr char startUtc[] = "2024-03-20T00:00:00";

struct Contact {
    int station{-1};
    fd::dynamics::CartesianState state;
};

CartesianPropagationConfig propagationConfig(const fd::simulation::EarthOrbitModel& model,
                                               double noiseDensity = 0) {
    CartesianPropagationConfig config;
    config.integrator.absoluteTolerance = 1e-11;
    config.integrator.relativeTolerance = 1e-11;
    config.integrator.initialStep = 30;
    config.integrator.maximumStep = 10;
    config.clock = model.clock;
    // NASA SNC: sqrt(PSD) is m/s^(3/2); the propagator uses km^2/s^3.
    config.accelerationDiffusion = std::pow(noiseDensity / 1000, 2) * Eigen::Matrix3d::Identity();
    return config;
}

Eigen::VectorXd priorSigmas(bool uncertainSrp, double positionSigma, double velocitySigma) {
    Eigen::VectorXd sigma(9);
    sigma << positionSigma / 1000, positionSigma / 1000, positionSigma / 1000,
             velocitySigma / 1000, velocitySigma / 1000, velocitySigma / 1000,
             1e-6, 1e-10, uncertainSrp ? srpSigma : 0;
    return sigma;
}

Eigen::VectorXd initialTruth(double mu) {
    Eigen::VectorXd state = Eigen::VectorXd::Zero(9);
    const double speed = std::sqrt(mu / 20000), inclination = 25 * std::numbers::pi / 180;
    state.head<6>() << 20000, 0, 0, 0, speed * std::cos(inclination), speed * std::sin(inclination);
    state[6] = 2e-7;
    state[7] = 2e-11;
    state[8] = 1;
    return state;
}

// Choose contacts from one nominal trajectory, so every run has the same contact schedule.
std::array<Contact, epochs> contactPlan(const fd::dynamics::SpiceEarthEnvironment& environment,
                                        double epoch, const CartesianPropagator& nominal) {
    std::vector<od::Station> stations;
    for (const char* name : {"DSS-43", "DSS-63", "DSS-14"})
        stations.push_back(
            od::buildStationFromKernel(name, od::stationNaifIdFromName(name), epoch));
    Eigen::VectorXd orbit = initialTruth(environment.earthMu()).head<6>();
    std::array<Contact, epochs> contacts;
    for (int k = 0; k < epochs; ++k) {
        const double time = k * step;
        if (k > 0)
            orbit = nominal(time - step, time, orbit).state;
        if ((time >= 7200 && time < 14400) || (time >= 21600 && time < 25200))
            continue;
        double best = std::sin(10 * std::numbers::pi / 180);
        for (std::size_t j = 0; j < stations.size(); ++j) {
            const auto state = environment.stationState(stations[j], epoch + time);
            const double elevation = (orbit.head<3>() - state.positionKm)
                                         .normalized()
                                         .dot(environment.stationUp(stations[j], epoch + time));
            if (elevation > best) {
                best = elevation;
                contacts[k] = {static_cast<int>(j), state};
            }
        }
    }
    return contacts;
}

void header(std::ostream& csv, int n) {
    csv << "run,time_s,tracking,station,position_error_m,position_sigma_max_m,nees,nees_orbit,"
           "nis,srp_true,srp_estimate,srp_sigma,observed_range_km,observed_rate_km_s";
    for (int j = 0; j < n; ++j)
        csv << ",error_" << j << ",sigma_" << j;
    csv << '\n' << std::setprecision(12);
}

void record(std::ostream& csv, std::size_t run, int k, const Contact& contact, const EKF& filter,
            const Eigen::VectorXd& truth, bool uncertainSrp, double nis,
            const Eigen::Vector2d& observed) {
    const int n = filter.state().size(), dof = n == 9 && uncertainSrp ? 9 : 8;
    const Eigen::VectorXd error = filter.state() - truth.head(n);
    const auto& covariance = filter.covariance();
    const double sigmaMax =
        std::sqrt(Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(covariance.topLeftCorner<3, 3>())
                      .eigenvalues()
                      .maxCoeff());
    // In the matched control alpha is known exactly: exclude that deterministic state from NEES.
    const double fullNees =
        fd::montecarlo::nees(error.head(dof), covariance.topLeftCorner(dof, dof));
    const double orbitNees =
        fd::montecarlo::nees(error.head<6>(), covariance.topLeftCorner<6, 6>());
    csv << run << ',' << k * step << ',' << (contact.station >= 0) << ',' << contact.station << ','
        << error.head<3>().norm() * 1000 << ',' << sigmaMax * 1000 << ',' << fullNees << ','
        << orbitNees << ',' << nis << ',' << truth[8] << ',' << (n == 9 ? filter.state()[8] : 1)
        << ',' << (n == 9 ? std::sqrt(covariance(8, 8)) : 0) << ',' << observed[0] << ','
        << observed[1];
    for (int j = 0; j < n; ++j)
        csv << ',' << error[j] << ',' << std::sqrt(covariance(j, j));
    csv << '\n';
}

void simulate(std::size_t run, std::uint64_t seed, bool uncertainSrp,
              const fd::dynamics::SpiceEarthEnvironment& environment, double epoch,
              const std::array<Contact, epochs>& contacts, std::ostream& eight, std::ostream& nine,
              double positionSigma, double velocitySigma, double noiseDensity) {
    const auto nominal = fd::simulation::nominalEarthOrbitModel();
    auto config = propagationConfig(nominal, noiseDensity);
    const auto force = fd::simulation::earthOrbitForces(environment, epoch, nominal);
    const CartesianPropagator propagate8(force, config);
    config.srpAcceleration = fd::simulation::earthOrbitSrp(environment, epoch, nominal);
    const CartesianPropagator propagate9(force, config);
    const Eigen::VectorXd sigma = priorSigmas(uncertainSrp, positionSigma, velocitySigma);
    const Eigen::MatrixXd prior = sigma.array().square().matrix().asDiagonal();
    std::mt19937_64 initialGenerator(seed), measurementGenerator(seed + 1);
    std::normal_distribution<double> initialNormal, measurementNormal;
    Eigen::VectorXd truth = initialTruth(environment.earthMu()), initial = truth;
    initial.head<8>() += fd::montecarlo::drawError(prior.topLeftCorner<8, 8>(), initialGenerator);
    truth[8] = uncertainSrp ? 1 + srpSigma * initialNormal(initialGenerator) : 1;
    if (truth[8] <= 0)
        throw std::runtime_error(
            "Gaussian SRP draw is nonphysical; no samples are silently clipped.");
    auto trueModel = nominal;
    trueModel.srpCr *= truth[8];
    const CartesianPropagator propagateTruth(
        fd::simulation::earthOrbitForces(environment, epoch, trueModel),
        propagationConfig(trueModel));
    const auto clock = fd::clocks::DSAC::fromParameters(nominal.clock);
    fd::clocks::ClockTruthSimulator clockTruth(clock, seed + 2);
    EKF filter8(StateLayout::OrbitClock), filter9(StateLayout::OrbitClockSrp);
    filter8.setInitialState(initial.head<8>(), prior.topLeftCorner<8, 8>(), 0);
    filter9.setInitialState(initial, prior, 0);
    const Eigen::Vector2d measurementSigma(.01, 1e-7);
    const Eigen::Matrix2d r = measurementSigma.array().square().matrix().asDiagonal();
    const double missing = std::numeric_limits<double>::quiet_NaN();
    for (int k = 0; k < epochs; ++k) {
        const double time = k * step;
        if (k > 0) {
            truth.head<6>() = propagateTruth(time - step, time, truth.head<6>()).state;
            const auto next = clockTruth.step({truth[6], truth[7]}, step);
            truth[6] = next.bias_s;
            truth[7] = next.fractional_frequency;
            filter8.predictTo(time, propagate8);
            filter9.predictTo(time, propagate9);
        }
        const double rangeNoise = measurementSigma[0] * measurementNormal(measurementGenerator);
        const double rateNoise = measurementSigma[1] * measurementNormal(measurementGenerator);
        Eigen::Vector2d observed = Eigen::Vector2d::Constant(missing);
        double nis8 = missing, nis9 = missing;
        const auto& contact = contacts[k];
        if (contact.station >= 0) {
            const auto model = [station = contact.state](double, const Eigen::VectorXd& x) {
                return fd::observations::geometricRadiometricPrediction(x, station);
            };
            observed = model(time, truth).value + Eigen::Vector2d(rangeNoise, rateNoise);
            nis8 = filter8.update(time, observed, r, model).normalizedInnovationSquared;
            nis9 = filter9.update(time, observed, r, model).normalizedInnovationSquared;
        }
        record(eight, run, k, contact, filter8, truth, uncertainSrp, nis8, observed);
        record(nine, run, k, contact, filter9, truth, uncertainSrp, nis9, observed);
    }
}

std::uint64_t integer(const std::string& text) {
    if (text.empty() || !std::all_of(text.begin(), text.end(),
                                     [](unsigned char c) { return c >= '0' && c <= '9'; }))
        throw std::invalid_argument("Expected an unsigned decimal integer.");
    return std::stoull(text);
}
} // namespace

int main(int argc, char** argv) {
    try {
        std::size_t runs = 100;
        std::uint64_t seed = 2026;
        double positionSigma = 100, velocitySigma = 1, noiseDensity = 1e-6;
        std::filesystem::path output = "montecarlo/data";
        std::filesystem::path kernels = std::filesystem::path(DEEPNAV_SOURCE_DIR) / "Kernels";
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (i + 1 >= argc)
                throw std::invalid_argument("Missing option value.");
            const std::string value = argv[++i];
            if (arg == "--runs")
                runs = integer(value);
            else if (arg == "--seed")
                seed = integer(value);
            else if (arg == "--position-sigma-m" || arg == "--velocity-sigma-m-s" ||
                     arg == "--acceleration-noise-density") {
                std::size_t consumed;
                const double parameter = std::stod(value, &consumed);
                if (consumed != value.size() || !std::isfinite(parameter) || parameter < 0 ||
                    (arg != "--acceleration-noise-density" && parameter == 0))
                    throw std::invalid_argument(
                        "Prior sigmas must be positive; noise density may be zero.");
                if (arg == "--position-sigma-m")
                    positionSigma = parameter;
                else if (arg == "--velocity-sigma-m-s")
                    velocitySigma = parameter;
                else
                    noiseDensity = parameter;
            } else if (arg == "--output")
                output = value;
            else if (arg == "--kernels")
                kernels = value;
            else
                throw std::invalid_argument(
                    "Usage: ekf_montecarlo [--runs N] [--seed N] [--position-sigma-m X] "
                    "[--velocity-sigma-m-s X] [--acceleration-noise-density X] "
                    "[--output DIR] [--kernels DIR]");
        }
        if (runs == 0 || runs > (std::numeric_limits<std::uint64_t>::max() - seed) / 3)
            throw std::invalid_argument("Run count must be positive and seeds must not overflow.");
        const fd::dynamics::SpiceEarthEnvironment environment(kernels);
        const double epoch = environment.epochTdb(startUtc);
        const auto nominal = fd::simulation::nominalEarthOrbitModel();
        const CartesianPropagator nominalPropagation(
            fd::simulation::earthOrbitForces(environment, epoch, nominal),
            propagationConfig(nominal));
        const auto contacts = contactPlan(environment, epoch, nominalPropagation);
        std::filesystem::create_directories(output);
        std::ofstream metadata(output / "study.json");
        metadata << std::setprecision(17) << "{\n  \"runs\": " << runs << ",\n  \"seed\": " << seed
                 << ",\n  \"start_utc\": \"" << startUtc << "\",\n  \"step_s\": " << step
                 << ",\n  \"duration_s\": " << duration << ",\n  \"srp_sigma\": " << srpSigma
                 << ",\n  \"acceleration_noise_density_m_s32\": " << noiseDensity
                 << ",\n  \"matched_nees_dof\": [8,8],\n  \"srp_nees_dof\": [8,9],\n"
                    "  \"initial_sigma\": [";
        const auto sigma = priorSigmas(false, positionSigma, velocitySigma);
        for (int j = 0; j < 8; ++j)
            metadata << (j ? "," : "") << sigma[j];
        metadata << "],\n  \"measurement_sigma\": [0.01,1e-7],\n"
                    "  \"stations\": [\"DSS-43\",\"DSS-63\",\"DSS-14\"],\n"
                    "  \"tracking\": \"fixed nominal visibility; gaps 2-4 h and 6-7 h\",\n"
                    "  \"force_model\": \"Earth + J2 + Sun/Moon + shadowed SRP\",\n"
                    "  \"process_noise_model\": \"isotropic SNC; integrated with orbital STM\",\n"
                    "  \"orbital_diffusion_km2_s3\": "
                 << std::pow(noiseDensity / 1000, 2)
                 << ",\n  \"srp_process_noise\": 0,\n"
                    "  \"clock_drift_per_s\": "
                 << nominal.clock.frequency_drift_per_s
                 << ",\n  \"clock_q_bias_s\": " << nominal.clock.q_bias_s
                 << ",\n  \"clock_q_frequency_per_s\": " << nominal.clock.q_frequency_per_s
                 << ",\n  \"seed_streams\": \"base+3*run: initial/SRP; +1: measurements; +2: "
                    "clock\"\n}\n";
        for (bool uncertain : {false, true}) {
            const std::string label = uncertain ? "srp" : "matched";
            std::ofstream eight(output / (label + "_8.csv")), nine(output / (label + "_9.csv"));
            if (!eight || !nine)
                throw std::runtime_error("Cannot create Monte Carlo CSVs.");
            header(eight, 8);
            header(nine, 9);
            for (std::size_t run = 0; run < runs; ++run) {
                simulate(run, seed + 3 * run, uncertain, environment, epoch, contacts, eight, nine,
                         positionSigma, velocitySigma, noiseDensity);
                if ((run + 1) % 10 == 0 || run + 1 == runs)
                    std::cout << label << ": " << run + 1 << '/' << runs << " paired runs\n"
                              << std::flush;
            }
            if (!eight || !nine)
                throw std::runtime_error("Monte Carlo CSV write failed.");
        }
        if (!metadata)
            throw std::runtime_error("Monte Carlo metadata write failed.");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
