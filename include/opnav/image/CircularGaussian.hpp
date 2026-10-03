#pragma once

#include "opnav/Types.hpp"
#include <cstdint>

namespace fd::opnav::image {

using Image = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using PixelMask = Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Parameters = Eigen::Matrix<double, 5, 1>;
using ParameterCovariance = Eigen::Matrix<double, 5, 5>;

struct CircularGaussian {
    PixelCoordinates center;
    double heightDn{1.0}; // Continuous peak intensity; total signal = 2*pi*h*sigma^2.
    double sigmaPixels{1.0}; // Standard deviation, not FWHM.
    double backgroundDn{0.0}; // Constant DN per pixel.
};

struct GaussianPixel {
    double dn;
    Parameters partials; // [sample, line, height, sigma, background].
};

// Integer sample/line coordinates denote pixel centers; boundaries are +/-0.5.
[[nodiscard]] GaussianPixel evaluateGaussianPixel(
    const CircularGaussian& model, PixelCoordinates pixel);
[[nodiscard]] Image renderCircularGaussian(
    ImageSize size, const CircularGaussian& model, PixelCoordinates origin = {});

struct DetectorNoise {
    double electronsPerDn{1.0};
    double readNoiseDn{0.0};
};

struct SyntheticImage {
    Image dn;
    Image varianceDn2; // Expected Poisson + read variance, before drawing noise.
};

[[nodiscard]] SyntheticImage simulateCircularGaussian(
    ImageSize size, const CircularGaussian& model, DetectorNoise noise,
    std::uint64_t seed, PixelCoordinates origin = {});

} // namespace fd::opnav::image
