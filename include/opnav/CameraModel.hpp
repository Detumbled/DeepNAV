#pragma once

#include "opnav/Types.hpp"
#include "opnav/DistortionModel.hpp"
#include <Eigen/Dense>
#include <memory>

namespace fd::opnav {

class CameraModel {
public:
    CameraModel(
        CameraIntrinsics intrinsics,
        std::shared_ptr<const DistortionModel> distortion);

    explicit CameraModel(CameraIntrinsics intrinsics)
        : CameraModel(std::move(intrinsics), std::make_shared<NoDistortion>()) {}

    [[nodiscard]]
    PixelCoordinates project(
        const Eigen::Vector3d& directionCameraToTarget) const;

    [[nodiscard]]
    const CameraIntrinsics& getIntrinsics() const noexcept { return intrinsics_; }

    [[nodiscard]]
    const DistortionModel& getDistortionModel() const noexcept { return *distortion_; }

private:
    // 1. Proiezione sul piano unitario Z = 1 (adimensionale)
    [[nodiscard]]
    NormalizedCoordinates projectIdeal(
        const Eigen::Vector3d& directionCameraToTarget) const;

    // 2. Passaggio da coordinate normalizzate distorte a pixel del sensore
    [[nodiscard]]
    PixelCoordinates normalizedToPixel(
        const NormalizedCoordinates& normalized) const;

    CameraIntrinsics intrinsics_;
    std::shared_ptr<const DistortionModel> distortion_;
};

} // namespace fd::opnav