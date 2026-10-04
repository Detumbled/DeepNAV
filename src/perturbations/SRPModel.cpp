#include "perturbations/SRP.hpp"
#include "perturbations/Eclipse.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::perturbations {
namespace {
void validateReflectivity(double cr) {
    if (!std::isfinite(cr) || cr < 0.0)
        throw std::invalid_argument(
            "SRP reflectivity coefficient must be finite and non-negative.");
}
void validatePositiveFinite(double value, const char* label) {
    if (!std::isfinite(value) || value <= 0.0)
        throw std::invalid_argument(std::string(label) + " must be positive and finite.");
}
} // namespace

SolarRadiationPressure::SolarRadiationPressure(std::string centralBody, std::string frame,
                                               double cr, double area_m2, double mass_kg)
    : centralBody_(std::move(centralBody)), frame_(std::move(frame)) {
    if (centralBody_.empty() || frame_.empty())
        throw std::invalid_argument("SRP central body and frame cannot be empty.");
    setReflectivity(cr);
    setArea(area_m2);
    setMass(mass_kg);
}

void SolarRadiationPressure::setReflectivity(double cr) {
    validateReflectivity(cr);
    cr_ = cr;
}
void SolarRadiationPressure::setArea(double area_m2) {
    validatePositiveFinite(area_m2, "SRP area");
    area_m2_ = area_m2;
}
void SolarRadiationPressure::setMass(double mass_kg) {
    validatePositiveFinite(mass_kg, "SRP mass");
    mass_kg_ = mass_kg;
}

AccelerationEvaluation
SolarRadiationPressure::evaluateAtSunPosition(const Eigen::Vector3d& scPosition,
                                              const Eigen::Vector3d& centralToSun) const {
    const Eigen::Vector3d relative = scPosition - centralToSun;
    const double radius = relative.norm();
    if (!scPosition.allFinite() || !centralToSun.allFinite() || !std::isfinite(radius) ||
        radius <= 0.0)
        throw std::invalid_argument(
            "SRP geometry must have finite nonzero Sun-spacecraft separation.");
    const Eigen::Vector3d unit = relative / radius;
    const double magnitude = kSolarPressure1AU_N_m2 * cr_ * (area_m2_ / mass_kg_) / 1000.0 *
                             std::pow(kAu_km / radius, 2.0);
    AccelerationEvaluation result{magnitude * unit,
                                  (magnitude / radius) *
                                      (Eigen::Matrix3d::Identity() - 3.0 * unit * unit.transpose()),
                                  Eigen::Matrix3d::Zero()};
    if (!result.acceleration.allFinite() || !result.positionJacobian.allFinite())
        throw std::runtime_error("SRP acceleration/partials overflowed.");
    return result;
}

AccelerationEvaluation SolarRadiationPressure::evaluateWithShadow(
    const Eigen::Vector3d& scPosition, const Eigen::Vector3d& centralToSun,
    double bodyRadiusKm, double sunRadiusKm) const {
    auto result = evaluateAtSunPosition(scPosition, centralToSun);
    const auto shadow = evaluateEclipse(scPosition, centralToSun, bodyRadiusKm, sunRadiusKm);
    result.positionJacobian = shadow.illumination * result.positionJacobian +
                              result.acceleration * shadow.positionGradient.transpose();
    result.acceleration *= shadow.illumination;
    return result;
}

} // namespace fd::perturbations
