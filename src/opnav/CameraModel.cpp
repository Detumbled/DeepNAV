#include "opnav/CameraModel.hpp"
#include <stdexcept>

namespace fd::opnav {

CameraModel::CameraModel(
    CameraIntrinsics intrinsics,
    std::shared_ptr<const DistortionModel> distortion)
    : intrinsics_(std::move(intrinsics)),
      distortion_(distortion ? std::move(distortion) : std::make_shared<NoDistortion>()) {}

NormalizedCoordinates CameraModel::projectIdeal(
    const Eigen::Vector3d& directionCameraToTarget) const {
    
    if (!directionCameraToTarget.allFinite()) {
        throw std::invalid_argument("Direction vector contains non-finite values.");
    }
    if (directionCameraToTarget.z() <= 0.0) {
        throw std::domain_error("Target direction must have positive Z component.");
    }

    return NormalizedCoordinates{
        directionCameraToTarget.x() / directionCameraToTarget.z(),
        directionCameraToTarget.y() / directionCameraToTarget.z()
    };
}

PixelCoordinates CameraModel::normalizedToPixel(
    const NormalizedCoordinates& normalized) const {

    const FocalPlaneCoordinates focal_plane_mm = intrinsics_.focalLength * normalized;

    // Conversione da millimetri a pixel tramite matrice di calibrazione
    const CameraIntrinsics::pixelOffset_t offset = 
        intrinsics_.focalPlaneToPixelMatrix_PixelPermm * focal_plane_mm;

    return PixelCoordinates{
        .sample = intrinsics_.principalPoint.sample + offset.x(),
        .line   = intrinsics_.principalPoint.line   + offset.y()
    };
}

PixelCoordinates CameraModel::project(
    const Eigen::Vector3d& directionCameraToTarget) const {

    const NormalizedCoordinates ideal = projectIdeal(directionCameraToTarget);
    const NormalizedCoordinates distorted = distortion_->distort(ideal);

    return normalizedToPixel(distorted);
}

} // namespace fd::opnav