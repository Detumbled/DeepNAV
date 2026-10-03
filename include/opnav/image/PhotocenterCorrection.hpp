#pragma once

#include "opnav/image/CentroidEstimator.hpp"

namespace fd::opnav::image {

struct PhotocenterOffset {
    Eigen::Vector2d pixels{Eigen::Vector2d::Zero()}; // Light center minus geometric center.
    Eigen::Matrix2d covariance{Eigen::Matrix2d::Zero()};
};

// Section 5.4, Eqs. 5.16-5.17: first moments of a supplied nonnegative,
// background-free projected brightness model, relative to its geometric center.
// This grid describes intrinsic body brightness, not the noisy observed PSF.
[[nodiscard]] PhotocenterOffset brightnessPhotocenterOffset(
    const Image& brightness, PixelCoordinates geometricCenter, PixelCoordinates origin = {});

// Assumes offset errors are independent of centroid errors. Returns a separate
// measurement so the raw photocenter stays available; apply exactly once.
[[nodiscard]] PixelMeasurement correctPhotocenter(
    const PixelMeasurement& photocenter, const PhotocenterOffset& offset);

} // namespace fd::opnav::image
