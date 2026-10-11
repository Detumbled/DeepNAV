#pragma once

#include "perturbations/Eclipse.hpp"
#include "perturbations/SRP.hpp"
#include <cmath>
#include <stdexcept>

namespace fd::perturbations {

// Opaque, one-sided Lambertian plate; no specular reflection. nu is the
// paper's coefficient (1-alpha)/3, not the diffuse energy fraction itself.
// At normal incidence this equals the cannonball law with CR=1+2*nu.
// Di Benedetto et al. (2009), Appendix Eqs. (6)-(8).
// The normal points outward from the illuminated face toward the Sun.
[[nodiscard]] inline AccelerationEvaluation flatPlateSRP(
    const Eigen::Vector3d& spacecraft, const Eigen::Vector3d& sun,
    const Eigen::Vector3d& normal, double areaM2, double massKg, double nu) {
    const Eigen::Vector3d q = spacecraft - sun;
    const double d = q.norm();
    if (!q.allFinite() || !spacecraft.allFinite() || !sun.allFinite() ||
        !std::isfinite(d) || d <= 0 || !normal.allFinite() ||
        std::abs(normal.norm()-1) > 1e-8 || !std::isfinite(areaM2) || areaM2 <= 0 ||
        !std::isfinite(massKg) || massKg <= 0 || !std::isfinite(nu) || nu < 0 || nu > 1.0/3)
        throw std::invalid_argument("Invalid flat-plate SRP geometry or optical parameters");
    AccelerationEvaluation result;
    const Eigen::Vector3d u = q/d;
    const double h = -normal.dot(u);
    if (h <= 0) return result; // Backside/bus forces require additional surfaces.
    const Eigen::Matrix3d U = (Eigen::Matrix3d::Identity()-u*u.transpose())/d;
    const Eigen::Vector3d g = -U*normal, momentum = u-2*nu*normal;
    const double k = SolarRadiationPressure::kSolarPressure1AU_N_m2 * areaM2/massKg/1000 *
        std::pow(SolarRadiationPressure::kAu_km/d,2);
    result.acceleration = k*h*momentum;
    result.positionJacobian = k*(h*U + momentum*(g-2*h*u/d).transpose());
    if (!result.acceleration.allFinite() || !result.positionJacobian.allFinite())
        throw std::runtime_error("Flat-plate SRP acceleration/partials overflowed");
    return result;
}

[[nodiscard]] inline AccelerationEvaluation flatPlateSRPWithShadow(
    const Eigen::Vector3d& spacecraft, const Eigen::Vector3d& sun,
    const Eigen::Vector3d& normal, double areaM2, double massKg, double nu,
    double bodyRadiusKm, double sunRadiusKm) {
    auto result = flatPlateSRP(spacecraft,sun,normal,areaM2,massKg,nu);
    const auto shadow = evaluateEclipse(spacecraft,sun,bodyRadiusKm,sunRadiusKm);
    result.positionJacobian = shadow.illumination*result.positionJacobian +
        result.acceleration*shadow.positionGradient.transpose();
    result.acceleration *= shadow.illumination;
    return result;
}

} // namespace fd::perturbations
