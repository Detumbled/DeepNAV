#include "opnav/core/GeometricLineOfSight.hpp"

#include <cmath>
#include <stdexcept>

namespace fd::opnav::core {

GeometricLineOfSight computeGeometricLineOfSight(
    const fd::dynamics::CartesianState& cameraState,
    const fd::dynamics::CartesianState& targetState) {

    if (!cameraState.allFinite() || !targetState.allFinite()) {
        throw std::invalid_argument(
            "Camera and target states must be finite.");
    }

    const Eigen::Vector3d vectorKm =
        targetState.positionKm - cameraState.positionKm;

    const double rangeKm = vectorKm.norm();

    if (!std::isfinite(rangeKm) || rangeKm <= 0.0) {
        throw std::domain_error(
            "Camera and target positions must be distinct.");
    }

    return GeometricLineOfSight {
        .vectorKm = vectorKm,
        .unitDirection = vectorKm / rangeKm,
        .rangeKm = rangeKm
    };
}

} // namespace fd::opnav::core
