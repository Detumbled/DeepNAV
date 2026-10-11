#pragma once

#include "perturbations/ForceModel.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::perturbations {

// Di Benedetto, Iess & Roth (2009), Eq. (4): body-fixed anisotropic
// emission, inverse mass scaling and Pu-238 decay. No solar-distance or
// eclipse factor. The signed coefficients are accelerations in km/s^2
// at referenceEpoch (seconds) and referenceMassKg.
// https://issfd.org/ISSFD_2009/OrbitDeterminationI/DiBENEDETTO.pdf
class ThermalRecoil {
public:
    ThermalRecoil(Eigen::Vector3d referenceBodyAcceleration, double referenceEpoch,
                  double referenceMassKg, double halfLifeSeconds = 87.74 * 365.25 * 86400)
        : body_(std::move(referenceBodyAcceleration)), epoch_(referenceEpoch),
          mass_(referenceMassKg), halfLife_(halfLifeSeconds) {
        if (!body_.allFinite() || !std::isfinite(epoch_) ||
            !std::isfinite(mass_) || mass_ <= 0 ||
            !std::isfinite(halfLife_) || halfLife_ <= 0)
            throw std::invalid_argument("Invalid thermal recoil reference parameters");
    }

    // Attitude and mass are prescribed independently of the propagated
    // Cartesian state, hence both state Jacobians are zero. Caller owns
    // the epoch convention, attitude coverage and mass history.
    [[nodiscard]] AccelerationEvaluation evaluate(
        double epoch, double massKg, const Eigen::Matrix3d& bodyToIntegration) const {
        if (!std::isfinite(epoch) || !std::isfinite(massKg) || massKg <= 0 ||
            !bodyToIntegration.allFinite() ||
            (bodyToIntegration.transpose() * bodyToIntegration -
             Eigen::Matrix3d::Identity()).norm() > 1e-8 ||
            std::abs(bodyToIntegration.determinant() - 1) > 1e-8)
            throw std::invalid_argument("Invalid thermal recoil epoch, mass or rotation");
        AccelerationEvaluation result;
        result.acceleration = bodyToIntegration * body_ *
            (mass_ / massKg * std::exp(-std::log(2.0) * ((epoch - epoch_) / halfLife_)));
        if (!result.acceleration.allFinite())
            throw std::runtime_error("Thermal recoil acceleration overflowed");
        return result;
    }
private:
    Eigen::Vector3d body_;
    double epoch_, mass_, halfLife_;
};

} // namespace fd::perturbations
