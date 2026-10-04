#pragma once

#include "perturbations/ForceModel.hpp"
#include <Eigen/Dense>

#include <string>
#include <vector>

namespace fd::perturbations {

[[nodiscard]] AccelerationFunction pointMassGravity(double muKm3PerSec2);
[[nodiscard]] AccelerationEvaluation thirdBodyGravity(double muKm3PerSec2,
                                                      const Eigen::Vector3d& spacecraftPosition,
                                                      const Eigen::Vector3d& centralToBody);

struct MassiveBody {
    std::string name;
    double mu{0.0};
};

class ThirdBodyGravity {
public:
    ThirdBodyGravity(std::string centralBody, std::string frame);

    void addBody(const MassiveBody& body);
    void clearBodies();

    [[nodiscard]] Eigen::Vector3d computeAcceleration(double tdb,
                                                      const Eigen::Vector3d& scPosition) const;
    [[nodiscard]] AccelerationEvaluation evaluate(double tdb,
                                                  const Eigen::Vector3d& scPosition) const;

private:
    std::string centralBody_;
    std::string frame_;
    std::vector<MassiveBody> perturbingBodies_;
};

} // namespace fd::perturbations
