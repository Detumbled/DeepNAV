#include "perturbations/Gravitational.hpp"
#include <cmath>
#include <stdexcept>

namespace fd::perturbations {
namespace {
void validateMu(double mu) {
    if (!std::isfinite(mu) || mu <= 0.0)
        throw std::invalid_argument("Gravitational parameter must be positive and finite.");
}
double separation(const Eigen::Vector3d& position) {
    const double radius = position.norm();
    if (!position.allFinite() || !std::isfinite(radius) || radius <= 0.0)
        throw std::invalid_argument("Gravity geometry must have finite nonzero separation.");
    return radius;
}
} // namespace

AccelerationFunction pointMassGravity(double mu) {
    validateMu(mu);
    return [mu](double, const Eigen::Vector3d& position, const Eigen::Vector3d&) {
        const double radius = separation(position);
        const Eigen::Vector3d unit = position / radius;
        const double factor = mu / radius / radius / radius;
        return AccelerationEvaluation{
            -factor * position,
            factor * (3.0 * unit * unit.transpose() - Eigen::Matrix3d::Identity()),
            Eigen::Matrix3d::Zero()};
    };
}

AccelerationEvaluation thirdBodyGravity(double mu, const Eigen::Vector3d& spacecraftPosition,
                                        const Eigen::Vector3d& centralToBody) {
    validateMu(mu);
    if (!spacecraftPosition.allFinite())
        throw std::invalid_argument("Spacecraft position must be finite.");
    const Eigen::Vector3d spacecraftToBody = centralToBody - spacecraftPosition;
    const double bodyRadius = separation(centralToBody), distance = separation(spacecraftToBody);
    const Eigen::Vector3d unit = spacecraftToBody / distance;
    const double factor = mu / distance / distance / distance;
    return {factor * spacecraftToBody -
                (mu / bodyRadius / bodyRadius) * (centralToBody / bodyRadius),
            factor * (3.0 * unit * unit.transpose() - Eigen::Matrix3d::Identity()),
            Eigen::Matrix3d::Zero()};
}

} // namespace fd::perturbations
