#include "perturbations/SRP.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <SpiceUsr.h>
#include <cmath>
#include <stdexcept>

namespace fd::perturbations {
namespace {
Eigen::Vector3d sunPositionRelativeToCentral(double tdb, const std::string& frame,
                                             const std::string& centralBody) {
    if (centralBody == "SUN" || centralBody == "10")
        return Eigen::Vector3d::Zero();
    od::SpiceErrorModeGuard actionGuard;
    SpiceDouble position[3] = {0.0, 0.0, 0.0}, lightTime = 0.0;
    spkpos_c("SUN", tdb, frame.c_str(), "NONE", centralBody.c_str(), position, &lightTime);
    od::throwIfSpiceFailed("Failed to compute Sun position for SRP.");
    return {position[0], position[1], position[2]};
}
} // namespace

Eigen::Vector3d
SolarRadiationPressure::computeAcceleration(double tdb, const Eigen::Vector3d& scPosition) const {
    return evaluate(tdb, scPosition).acceleration;
}

AccelerationEvaluation SolarRadiationPressure::evaluate(double tdb,
                                                        const Eigen::Vector3d& scPosition) const {
    if (!std::isfinite(tdb) || !scPosition.allFinite())
        throw std::invalid_argument("SRP acceleration inputs must be finite.");
    return evaluateAtSunPosition(scPosition,
                                 sunPositionRelativeToCentral(tdb, frame_, centralBody_));
}

} // namespace fd::perturbations
