#include "Clocks/ClockModel.hpp"
#include "dynamics/SpiceEarthEnvironment.hpp"
#include "filters/CartesianPropagator.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/J2.hpp"
#include "perturbations/SRP.hpp"
#include "simulation/EarthOrbitScenario.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace

int main() {
    try {
        using namespace fd::simulation;
        const fd::dynamics::SpiceEarthEnvironment environment(
            std::filesystem::path(DEEPNAV_SOURCE_DIR) / "Kernels");
        const double epoch = environment.epochTdb("2024-03-20T00:00:00");
        const auto truth = nominalEarthOrbitModel();
        const auto mismatch = namedModelMismatch("combined");
        const auto estimate = estimatorModel(truth, mismatch);
        const auto trueForce = earthOrbitForces(environment, epoch, truth);
        const auto estimatedForce = earthOrbitForces(environment, epoch, estimate);
        const Eigen::Vector3d position(20000, 4000, 3000), velocity(0, 4, 1);
        const auto actual = trueForce(0, position, velocity);
        const auto sun = environment.sunPosition(epoch), moon = environment.moonPosition(epoch);
        const auto gravity = fd::perturbations::pointMassGravity(environment.earthMu());
        const auto oblate =
            fd::perturbations::j2Gravity(environment.earthMu(), environment.earthRadius(),
                                         kNominalEarthJ2, position, environment.earthPole(epoch));
        const auto solar = fd::perturbations::thirdBodyGravity(environment.sunMu(), position, sun);
        const fd::perturbations::SolarRadiationPressure srp("EARTH", "J2000", truth.srpCr,
                                                            truth.srpAreaM2, truth.srpMassKg);
        const auto radiation = srp.evaluateWithShadow(position, sun, environment.earthRadius(),
                                                      environment.sunRadius());
        const Eigen::Vector3d delta = moon - position;
        const Eigen::Vector3d referenceMoon =
            environment.moonMu() *
            (delta / std::pow(delta.norm(), 3) - moon / std::pow(moon.norm(), 3));
        const Eigen::Vector3d recoveredMoon =
            actual.acceleration - gravity(0, position, velocity).acceleration -
            oblate.acceleration - solar.acceleration - radiation.acceleration;
        require((recoveredMoon - referenceMoon).norm() / referenceMoon.norm() < 1e-9,
                "Scenario lunar force is missing or has an incorrect origin.");
        const Eigen::Vector3d forceError =
            estimatedForce(0, position, velocity).acceleration - actual.acceleration;
        const Eigen::Vector3d expectedError =
            (mismatch.srpCrScale * mismatch.srpAreaMassScale - 1) * radiation.acceleration;
        require((forceError - expectedError).norm() / expectedError.norm() < 1e-7,
                "SRP mismatch leaked into other scenario forces.");
        const auto matched = earthOrbitForces(environment, epoch, estimatorModel(truth, {}));
        require((matched(0, position, velocity).acceleration - actual.acceleration).norm() == 0,
                "Matched model does not reproduce truth forces.");

        fd::filters::CartesianPropagationConfig truthConfig, estimateConfig;
        truthConfig.clock = truth.clock;
        estimateConfig.clock = estimate.clock;
        const fd::filters::CartesianPropagator propagateTruth(gravity, truthConfig),
            propagateEstimate(gravity, estimateConfig);
        Eigen::VectorXd initial(8);
        initial << position, velocity, 2e-7, 2e-11;
        constexpr double duration = 12000;
        const auto a = propagateTruth(0, duration, initial),
                   b = propagateEstimate(0, duration, initial);
        require((a.state.head<6>() - b.state.head<6>()).norm() == 0,
                "Clock mismatch changed orbital truth dynamics.");
        const double driftOffset = mismatch.clockDriftOffsetPerDay / 86400;
        require(std::abs((b.state[6] - a.state[6]) / (.5 * driftOffset * duration * duration) - 1) <
                        1e-10 &&
                    std::abs((b.state[7] - a.state[7]) / (driftOffset * duration) - 1) < 1e-10,
                "Clock drift mismatch does not produce analytic holdover error.");
        require(std::abs(b.processCovariance(6, 6) / a.processCovariance(6, 6) -
                         mismatch.clockNoiseScale) < 1e-12,
                "Clock noise mismatch did not scale estimator process covariance.");
        bool rejected = false;
        try {
            (void)estimatorModel(truth, {1, -1, 0, 1});
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "Invalid mismatch must fail before propagation.");
        std::cout << "Lunar force composition, fixed truth, SRP and clock mismatch holdover checks "
                     "passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
