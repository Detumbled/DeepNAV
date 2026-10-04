#include "perturbations/J2.hpp"
#include <cmath>
#include <stdexcept>

namespace fd::perturbations {
AccelerationEvaluation j2Gravity(double mu, double equatorialRadius, double j2,
                                 const Eigen::Vector3d& position, const Eigen::Vector3d& pole) {
    const double radius = position.norm();
    if (!std::isfinite(mu) || mu <= 0 || !std::isfinite(equatorialRadius) ||
        equatorialRadius <= 0 || !std::isfinite(j2) || !position.allFinite() ||
        !std::isfinite(radius) || radius <= 0 || !pole.allFinite() ||
        std::abs(pole.squaredNorm() - 1) > 1e-12)
        throw std::invalid_argument(
            "J2 requires positive mu/radius, finite geometry and a unit pole.");
    const double radius2 = radius * radius, z = position.dot(pole);
    const double shape = 5 * z * z / radius2 - 1;
    const double factor =
        1.5 * j2 * (mu / radius2) * std::pow(equatorialRadius / radius, 2) / radius;
    const Eigen::Vector3d direction = shape * position - 2 * z * pole;
    const Eigen::Vector3d shapeGradient =
        10 * z / radius2 * pole - 10 * z * z / (radius2 * radius2) * position;
    AccelerationEvaluation result;
    result.acceleration = factor * direction;
    result.positionJacobian =
        factor * (shape * Eigen::Matrix3d::Identity() + position * shapeGradient.transpose() -
                  2 * pole * pole.transpose() - 5 / radius2 * direction * position.transpose());
    if (!result.acceleration.allFinite() || !result.positionJacobian.allFinite())
        throw std::runtime_error("J2 acceleration/partials overflowed.");
    return result;
}
} // namespace fd::perturbations
