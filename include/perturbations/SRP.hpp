#pragma once

#include "perturbations/ForceModel.hpp"
#include <Eigen/Dense>

#include <string>

namespace fd::perturbations {

class SolarRadiationPressure {
public:
    SolarRadiationPressure(std::string centralBody, std::string frame, double cr, double area_m2,
                           double mass_kg);

    void setReflectivity(double cr);
    void setArea(double area_m2);
    void setMass(double mass_kg);

    [[nodiscard]] Eigen::Vector3d computeAcceleration(double tdb,
                                                      const Eigen::Vector3d& scPosition) const;
    [[nodiscard]] AccelerationEvaluation evaluate(double tdb,
                                                  const Eigen::Vector3d& scPosition) const;
    // Geometry-only evaluation reuses the same force law with caller-supplied ephemerides.
    [[nodiscard]] AccelerationEvaluation
    evaluateAtSunPosition(const Eigen::Vector3d& scPosition,
                          const Eigen::Vector3d& centralToSun) const;
    [[nodiscard]] AccelerationEvaluation
    evaluateWithShadow(const Eigen::Vector3d& scPosition, const Eigen::Vector3d& centralToSun,
                       double bodyRadiusKm, double sunRadiusKm) const;

    static constexpr double kSolarPressure1AU_N_m2 = 4.5605e-6;
    static constexpr double kAu_km = 149597870.7;

private:
    std::string centralBody_;
    std::string frame_;
    double cr_{1.0};
    double area_m2_{1.0};
    double mass_kg_{1.0};
};

} // namespace fd::perturbations
