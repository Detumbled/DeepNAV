#include "perturbations/Eclipse.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace fd::perturbations {
EclipseEvaluation evaluateEclipse(const Eigen::Vector3d& spacecraft, const Eigen::Vector3d& sun,
                                  double bodyRadiusKm, double sunRadiusKm) {
    const Eigen::Vector3d toBody = -spacecraft, toSun = sun - spacecraft;
    const double bodyDistance = toBody.norm(), sunDistance = toSun.norm();
    if (!spacecraft.allFinite() || !sun.allFinite() || !std::isfinite(bodyRadiusKm) ||
        !std::isfinite(sunRadiusKm) || !std::isfinite(bodyDistance) ||
        !std::isfinite(sunDistance) || bodyRadiusKm <= 0 || sunRadiusKm <= 0 ||
        bodyDistance <= bodyRadiusKm || sunDistance <= sunRadiusKm)
        throw std::invalid_argument("Eclipse geometry requires finite radii and an "
                                    "observer outside both bodies.");
    const Eigen::Vector3d u = toBody / bodyDistance, v = toSun / sunDistance;
    const double sunRatio = sunRadiusKm / sunDistance, bodyRatio = bodyRadiusKm / bodyDistance;
    const double a = std::asin(sunRatio), b = std::asin(bodyRatio);
    const double separationSine = u.cross(v).norm();
    const double d = std::atan2(separationSine, u.dot(v));
    EclipseEvaluation result;
    if (d >= a + b)
        return result;
    if (b >= d + a) {
        result.illumination = 0;
        result.flag = EclipseFlag::Umbra;
        return result;
    }
    result.flag = EclipseFlag::Penumbra;
    const Eigen::Vector3d gradA = sunRatio / (sunDistance * std::sqrt(1 - sunRatio * sunRatio)) * v;
    const Eigen::Vector3d gradB =
        bodyRatio / (bodyDistance * std::sqrt(1 - bodyRatio * bodyRatio)) * u;
    if (a >= d + b) {
        // Annular geometry: the occulting disk is entirely inside the solar disk.
        result.illumination = std::clamp(1 - b * b / (a * a), 0.0, 1.0);
        result.positionGradient = -2 * b / (a * a) * gradB + 2 * b * b / (a * a * a) * gradA;
        return result;
    }
    const double x = std::clamp((d * d + a * a - b * b) / (2 * d * a), -1.0, 1.0);
    const double y = std::clamp((d * d + b * b - a * a) / (2 * d * b), -1.0, 1.0);
    const double root =
        std::sqrt(std::max(0.0, (-d + a + b) * (d + a - b) * (d - a + b) * (d + a + b)));
    const double overlap = a * a * std::acos(x) + b * b * std::acos(y) - .5 * root;
    const double solarArea = std::numbers::pi * a * a;
    result.illumination = std::clamp(1 - overlap / solarArea, 0.0, 1.0);
    const Eigen::Vector3d gradD =
        ((v - u.dot(v) * u) / bodyDistance + (u - u.dot(v) * v) / sunDistance) / separationSine;
    // Disk-overlap derivatives: radius derivatives are arc lengths; separation
    // derivative is -chord.
    const Eigen::Vector3d gradOverlap =
        2 * a * std::acos(x) * gradA + 2 * b * std::acos(y) * gradB - root / d * gradD;
    result.positionGradient = -gradOverlap / solarArea + 2 * overlap / (solarArea * a) * gradA;
    return result;
}
} // namespace fd::perturbations
