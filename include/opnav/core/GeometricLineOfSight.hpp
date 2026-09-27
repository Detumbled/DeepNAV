#pragma once

#include "dynamics/CartesianState.hpp"

#include <Eigen/Core>

namespace fd::opnav::core {

struct GeometricLineOfSight {
    Eigen::Vector3d vectorKm;
    Eigen::Vector3d unitDirection;
    double rangeKm {0.0};
};

[[nodiscard]]
GeometricLineOfSight computeGeometricLineOfSight(
    const fd::dynamics::CartesianState& cameraState,
    const fd::dynamics::CartesianState& targetState);

} // namespace fd::opnav::core
