#pragma once

#include <Eigen/Geometry>

namespace fd::perturbations {

enum class EclipseFlag { Sunlit = 0, Penumbra = 1, Umbra = 2 };

struct EclipseEvaluation {
    double illumination{1.0};
    EclipseFlag flag{EclipseFlag::Sunlit};
    Eigen::Vector3d positionGradient{Eigen::Vector3d::Zero()};
};

// Spherical bodies, uniform solar disk; overlap is evaluated in angular disk
// coordinates. All positions are relative to the occulting body's center, in km
// in one common frame.
[[nodiscard]] EclipseEvaluation evaluateEclipse(const Eigen::Vector3d& spacecraft,
                                                const Eigen::Vector3d& sun, double bodyRadiusKm,
                                                double sunRadiusKm);

} // namespace fd::perturbations
