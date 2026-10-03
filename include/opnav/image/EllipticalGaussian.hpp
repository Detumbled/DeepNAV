#pragma once

#include "opnav/image/CentroidEstimator.hpp"

namespace fd::opnav::image {

using EllipticalParameters = Eigen::Matrix<double, 7, 1>;
using EllipticalParameterCovariance = Eigen::Matrix<double, 7, 7>;

struct EllipticalGaussian {
    PixelCoordinates center;
    double heightDn{1};
    double sigmaMajorPixels{1};
    double sigmaMinorPixels{1};
    // Radians from +sample toward +line (clockwise on an image displayed line-down).
    double angleRadians{0};
    double backgroundDn{0};
};

struct EllipticalGaussianPixel {
    double dn;
    // [sample, line, height, sigmaMajor, sigmaMinor, angle, background].
    EllipticalParameters partials;
};

// Integrates over the actual square pixel using 16x16 Gauss-Legendre quadrature.
// Widths must be >= 0.2 pixels; integer coordinates are pixel centers.
[[nodiscard]] EllipticalGaussianPixel evaluateEllipticalGaussianPixel(
    const EllipticalGaussian& model, PixelCoordinates pixel);
[[nodiscard]] Image renderEllipticalGaussian(
    ImageSize size, const EllipticalGaussian& model, PixelCoordinates origin = {});

struct EllipticalGaussianFitResult {
    FitStatus status{FitStatus::NonConverged};
    EllipticalGaussian model; // Canonical major >= minor; angle in [-pi/2, pi/2).
    std::optional<PixelMeasurement> measurement;
    // Physical parameter order as above. Absent for a nearly circular model,
    // whose angle is undefined; the centroid covariance remains available.
    std::optional<EllipticalParameterCovariance> parameterCovariance;
    double chiSquared{0};
    std::size_t degreesOfFreedom{0};
    std::size_t iterations{0};
};

// Same fixed-variance WLS, masks, origins, status, and options as the circular fit.
// Positive-definite shape is optimized without an angle singularity at circularity.
[[nodiscard]] EllipticalGaussianFitResult fitEllipticalGaussian(
    const Image& dn, const Image& varianceDn2, PixelCoordinates origin = {},
    const PixelMask* mask = nullptr,
    std::optional<EllipticalGaussian> initial = std::nullopt,
    GaussianFitOptions options = {});

[[nodiscard]] EllipticalGaussianFitResult fitBrightestEllipticalGaussian(
    const Image& dn, const Image& varianceDn2, int windowRadius = 5,
    PixelCoordinates origin = {}, const PixelMask* mask = nullptr,
    GaussianFitOptions options = {});

} // namespace fd::opnav::image
