#include "simulation/EarthOrbitScenario.hpp"
#include "Clocks/DSAC.hpp"
#include "dynamics/SpiceEarthEnvironment.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/J2.hpp"
#include "perturbations/SRP.hpp"
#include <cmath>
#include <stdexcept>

namespace fd::simulation {
EarthOrbitModel nominalEarthOrbitModel() {
    EarthOrbitModel model;
    model.clock = fd::clocks::DSAC::shortTermWhiteFmBaseline().parameters();
    return model;
}

ModelMismatch namedModelMismatch(const std::string& name) {
    if (name == "matched")
        return {};
    if (name == "srp")
        return {.95, .90, 0, 1};
    if (name == "clock")
        return {1, 1, 3e-12, .5};
    if (name == "combined")
        return {.95, .90, 3e-12, .5};
    throw std::invalid_argument("--mismatch must be matched, srp, clock or combined.");
}

EarthOrbitModel estimatorModel(const EarthOrbitModel& truth, const ModelMismatch& mismatch) {
    if (!std::isfinite(mismatch.srpCrScale) || mismatch.srpCrScale < 0 ||
        !std::isfinite(mismatch.srpAreaMassScale) || mismatch.srpAreaMassScale <= 0 ||
        !std::isfinite(mismatch.clockDriftOffsetPerDay) ||
        !std::isfinite(mismatch.clockNoiseScale) || mismatch.clockNoiseScale < 0)
        throw std::invalid_argument(
            "Mismatch scales must be finite and nonnegative (area/mass positive).");
    EarthOrbitModel model = truth;
    model.srpCr *= mismatch.srpCrScale;
    model.srpAreaM2 *= mismatch.srpAreaMassScale;
    model.clock.frequency_drift_per_s += mismatch.clockDriftOffsetPerDay / 86400;
    model.clock.q_bias_s *= mismatch.clockNoiseScale;
    model.clock.q_frequency_per_s *= mismatch.clockNoiseScale;
    // Existing modules validate physical parameters and detect overflow before simulation starts.
    (void)fd::perturbations::SolarRadiationPressure("EARTH", "J2000", model.srpCr, model.srpAreaM2,
                                                    model.srpMassKg);
    (void)fd::clocks::ClockModel(model.clock);
    return model;
}

fd::perturbations::AccelerationFunction
earthOrbitForces(const fd::dynamics::SpiceEarthEnvironment& environment, double startEpoch,
                 const EarthOrbitModel& model, bool twoBodyOnly) {
    using namespace fd::perturbations;
    auto gravity = pointMassGravity(environment.earthMu());
    if (twoBodyOnly)
        return gravity;
    const SolarRadiationPressure srp("EARTH", "J2000", model.srpCr, model.srpAreaM2,
                                     model.srpMassKg);
    // environment must outlive the returned callback. Forces share Earth origin and J2000 axes.
    return [&environment, startEpoch, gravity, srp](double elapsed, const Eigen::Vector3d& position,
                                                    const Eigen::Vector3d& velocity) {
        const double epoch = startEpoch + elapsed;
        const auto sun = environment.sunPosition(epoch);
        auto force = gravity(elapsed, position, velocity);
        const auto oblate = j2Gravity(environment.earthMu(), environment.earthRadius(),
                                      kNominalEarthJ2, position, environment.earthPole(epoch));
        const auto solar = thirdBodyGravity(environment.sunMu(), position, sun);
        const auto lunar =
            thirdBodyGravity(environment.moonMu(), position, environment.moonPosition(epoch));
        const auto radiation = srp.evaluateWithShadow(position, sun, environment.earthRadius(),
                                                      environment.sunRadius());
        force.acceleration +=
            oblate.acceleration + solar.acceleration + lunar.acceleration + radiation.acceleration;
        force.positionJacobian += oblate.positionJacobian + solar.positionJacobian +
                                  lunar.positionJacobian + radiation.positionJacobian;
        return force;
    };
}
} // namespace fd::simulation
