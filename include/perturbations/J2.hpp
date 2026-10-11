#pragma once

#include "perturbations/ForceModel.hpp"

namespace fd::perturbations {

// J2 correction only: add to point-mass gravity. The unit pole and position share a frame.
// mu: km^3/s^2, equatorial radius/position: km, J2: unnormalized dimensionless coefficient.
[[nodiscard]] AccelerationEvaluation
j2Gravity(double mu, double equatorialRadius, double j2, const Eigen::Vector3d& position,
          const Eigen::Vector3d& pole = Eigen::Vector3d::UnitZ());

// Minus the central body's reaction acceleration due to its J2 interaction
// with a prescribed moon (Moyer 4.4.6). Independent of spacecraft state.
[[nodiscard]] AccelerationEvaluation
j2OriginCorrection(double centralGM, double moonGM, double radius, double j2,
                   const Eigen::Vector3d& moonPosition, const Eigen::Vector3d& pole);

} // namespace fd::perturbations
