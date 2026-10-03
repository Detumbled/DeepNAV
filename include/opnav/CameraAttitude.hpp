#pragma once

#include <Eigen/Core>

namespace fd::opnav {

// Passive rotation: directionCamera = C * directionInertial.
// Supply the camera attitude at the midpoint of the exposure.
class CameraAttitude {
public:
    explicit CameraAttitude(const Eigen::Matrix3d& inertialToCamera);

    [[nodiscard]] Eigen::Vector3d toCamera(const Eigen::Vector3d& directionInertial) const;
    [[nodiscard]] const Eigen::Matrix3d& matrix() const noexcept { return inertialToCamera_; }

private:
    Eigen::Matrix3d inertialToCamera_;
};

} // namespace fd::opnav
