#include "opnav/CameraAttitude.hpp"

#include <Eigen/LU>
#include <cmath>
#include <stdexcept>

namespace fd::opnav {

CameraAttitude::CameraAttitude(const Eigen::Matrix3d& inertialToCamera)
    : inertialToCamera_(inertialToCamera) {
    // Validate once so the same attitude can be reused for every target.
    if (!inertialToCamera_.allFinite()
        || (inertialToCamera_ * inertialToCamera_.transpose()
            - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() > 1e-12
        || std::abs(inertialToCamera_.determinant() - 1.0) > 1e-12) {
        throw std::invalid_argument("Camera attitude must be a finite proper rotation.");
    }
}

Eigen::Vector3d CameraAttitude::toCamera(const Eigen::Vector3d& directionInertial) const {
    if (!directionInertial.allFinite() || directionInertial.isZero(0.0)) {
        throw std::invalid_argument("Inertial direction must be finite and nonzero.");
    }
    const Eigen::Vector3d result = inertialToCamera_ * directionInertial;
    if (!result.allFinite()) throw std::overflow_error("Camera direction overflowed.");
    return result;
}

} // namespace fd::opnav
