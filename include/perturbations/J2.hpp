#pragma once

#include "perturbations/ForceModel.hpp"

namespace fd::perturbations {

// J2 correction only: add to point-mass gravity. The unit pole and position share a frame.
// mu: km^3/s^2, equatorial radius/position: km, J2: unnormalized dimensionless coefficient.
[[nodiscard]] AccelerationEvaluation
j2Gravity(double mu, double equatorialRadius, double j2, const Eigen::Vector3d& position,
          const Eigen::Vector3d& pole = Eigen::Vector3d::UnitZ());

} // namespace fd::perturbations
