#pragma once

#include "opnav/image/CircularGaussian.hpp"
#include <cstddef>
#include <optional>

namespace fd::opnav::image {

struct PixelMeasurement {
    PixelCoordinates center;
    Eigen::Matrix2d covariance; // sample/line covariance in pixels^2.
};

enum class FitStatus { Converged, NoSignal, NonConverged, Singular };

struct GaussianFitOptions {
    std::size_t maxIterations{80};
    double stepTolerance{1e-7};
    double minimumHeightSnr{5.0};
};

struct GaussianFitResult {
    FitStatus status{FitStatus::NonConverged};
    CircularGaussian model;
    std::optional<PixelMeasurement> measurement; // Only populated on success.
    std::optional<ParameterCovariance> parameterCovariance;
    double chiSquared{0};
    std::size_t degreesOfFreedom{0};
    std::size_t iterations{0};
};

// One isolated source in a calibrated local window; variance is known DN^2.
// A nonzero mask entry excludes that pixel (including saturation/bad pixels).
// Origin gives the full-image coordinate of window pixel (0,0).
// Covariance is local linearized WLS covariance with fixed supplied variances.
[[nodiscard]] GaussianFitResult fitCircularGaussian(
    const Image& dn, const Image& varianceDn2, PixelCoordinates origin = {},
    const PixelMask* mask = nullptr,
    std::optional<CircularGaussian> initial = std::nullopt,
    GaussianFitOptions options = {});

} // namespace fd::opnav::image
