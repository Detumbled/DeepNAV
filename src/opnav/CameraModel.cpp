#include "opnav/CameraModel.hpp"
#include <stdexcept>
#include <Eigen/LU>

namespace fd::opnav {

Eigen::Vector3d CameraModel::unproject(PixelCoordinates pixel) const {
    if(!std::isfinite(pixel.sample) || !std::isfinite(pixel.line)) throw std::invalid_argument("Invalid detector pixel.");
    const Eigen::Matrix2d focal=intrinsics_.focalLength*intrinsics_.focalPlaneToPixelMatrix_PixelPermm;
    if(!focal.allFinite() || std::abs(focal.determinant())<1e-15) throw std::invalid_argument("Singular camera calibration.");
    const Eigen::Vector2d observed=focal.inverse()*Eigen::Vector2d(pixel.sample-intrinsics_.principalPoint.sample,pixel.line-intrinsics_.principalPoint.line);
    Eigen::Vector2d ideal=observed;
    if(!dynamic_cast<const NoDistortion*>(distortion_.get())) {
        bool converged=false;
        for(int iteration=0;iteration<20;++iteration) {
            const Eigen::Vector2d residual=distortion_->distort(ideal)-observed;
            if(residual.norm()<1e-14){converged=true;break;}
            Eigen::Matrix2d jac;
            for(int j=0;j<2;++j) {Eigen::Vector2d step=Eigen::Vector2d::Zero();step[j]=1e-7;
                jac.col(j)=(distortion_->distort(ideal+step)-distortion_->distort(ideal-step))/(2e-7);}
            if(!jac.allFinite() || std::abs(jac.determinant())<1e-15)break;
            ideal-=jac.fullPivLu().solve(residual);
        }
        if(!converged)throw std::runtime_error("Inverse camera distortion did not converge.");
    }
    return Eigen::Vector3d(ideal.x(),ideal.y(),1).normalized();
}

PixelCoordinates CameraModel::project(
    const Eigen::Vector3d& directionInertial,
    const CameraAttitude& attitude) const {
    return project(attitude.toCamera(directionInertial));
}

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
