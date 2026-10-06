#include "observations/radiometric/PropagationCorrections.hpp"
#include "perturbations/Shapiro.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace fd::observations::radiometric {
double pointMassShapiroSeconds(const LegSolution& leg, const dynamics::StateProvider& body,
                               double gm) {
    if (!std::isfinite(gm) || gm <= 0)
        throw std::invalid_argument("Shapiro GM must be positive and finite.");
    const auto b1 = body.stateAt(opnav::TdbEpoch{leg.emissionTdb});
    const auto b2 = body.stateAt(opnav::TdbEpoch{leg.receptionTdb});
    const double r1 = (leg.emitter.positionKm - b1.positionKm).norm();
    const double r2 = (leg.receiver.positionKm - b2.positionKm).norm();
    const double distance = (leg.receiver.positionKm - leg.emitter.positionKm).norm();
    const double denominator = r1 + r2 - distance;
    if (!b1.allFinite() || !b2.allFinite() || !std::isfinite(denominator) || denominator <= 0)
        throw std::runtime_error("Singular Shapiro geometry.");
    // First-order point-mass GR delay (gamma=1). Not a conjunction-grade model.
    return fd::perturbations::computeShapiroTimeDelay(r1, r2, distance, gm);
}

double zenithPathDelaySeconds(double zenithPathM, double elevation) {
    if (!std::isfinite(zenithPathM) || zenithPathM < 0 || !std::isfinite(elevation) ||
        elevation <= 0 || elevation > std::numbers::pi / 2)
        throw std::invalid_argument("Zenith path must be nonnegative, elevation in (0,pi/2].");
    return zenithPathM / std::sin(elevation) / (1000 * speedOfLightKmPerSecond);
}
} // namespace fd::observations::radiometric
