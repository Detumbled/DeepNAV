#pragma once

#include <Eigen/Core>
#include <functional>
#include <vector>

namespace fd::perturbations {

// Acceleration in km/s^2; partials with respect to Cartesian km and km/s.
struct AccelerationEvaluation {
    Eigen::Vector3d acceleration{Eigen::Vector3d::Zero()};
    Eigen::Matrix3d positionJacobian{Eigen::Matrix3d::Zero()};
    Eigen::Matrix3d velocityJacobian{Eigen::Matrix3d::Zero()};
};
using AccelerationFunction = std::function<AccelerationEvaluation(
    double epoch, const Eigen::Vector3d& position, const Eigen::Vector3d& velocity)>;

// Sum forces and their partials using one consistent frame, origin and epoch.
[[nodiscard]] AccelerationFunction sumAccelerations(std::vector<AccelerationFunction> forces);

} // namespace fd::perturbations
