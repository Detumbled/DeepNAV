#pragma once

#include "opnav/image/CircularGaussian.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace fd::opnav::image::detail {

inline void validate(const CircularGaussian& m) {
    if (!std::isfinite(m.center.sample) || !std::isfinite(m.center.line)
        || !std::isfinite(m.heightDn) || m.heightDn <= 0
        || !std::isfinite(m.sigmaPixels) || m.sigmaPixels <= 0
        || !std::isfinite(m.backgroundDn))
        throw std::invalid_argument("Gaussian parameters must be finite, with positive height and width.");
}

inline void validateOrigin(PixelCoordinates origin) {
    if (!std::isfinite(origin.sample) || !std::isfinite(origin.line))
        throw std::invalid_argument("Image origin must be finite.");
}

struct AxisIntegral { double value, centerPartial, widthPartial; };

inline AxisIntegral integrateAxis(double pixel, double center, double sigma) {
    const double lo = pixel - 0.5 - center, hi = pixel + 0.5 - center;
    const double a = lo / sigma, b = hi / sigma;
    constexpr double root2 = std::numbers::sqrt2;
    // erfc avoids cancellation in same-sign tails.
    const double probability = a >= 0 ? std::erfc(a/root2) - std::erfc(b/root2)
        : b <= 0 ? std::erfc(-b/root2) - std::erfc(-a/root2)
        : std::erf(b/root2) - std::erf(a/root2);
    const double value = sigma * std::sqrt(std::numbers::pi / 2) * probability;
    const double ea = std::exp(-0.5*a*a), eb = std::exp(-0.5*b*b);
    // Avoid infinity * zero for distant pixels or very narrow models.
    const double aea = ea == 0 ? 0 : a*ea, beb = eb == 0 ? 0 : b*eb;
    return {value, ea - eb, value / sigma + aea - beb};
}

inline GaussianPixel combine(const CircularGaussian& m, AxisIntegral x, AxisIntegral y) {
    const double shape = x.value * y.value;
    GaussianPixel result;
    result.dn = m.backgroundDn + m.heightDn * shape;
    result.partials << m.heightDn*x.centerPartial*y.value,
        m.heightDn*x.value*y.centerPartial, shape,
        m.heightDn*(x.widthPartial*y.value + x.value*y.widthPartial), 1.0;
    if (!std::isfinite(result.dn) || !result.partials.allFinite())
        throw std::overflow_error("Gaussian pixel calculation overflowed.");
    return result;
}

inline void fillAxes(const CircularGaussian& m, Eigen::Index rows, Eigen::Index cols,
                    PixelCoordinates origin, std::vector<AxisIntegral>& x,
                    std::vector<AxisIntegral>& y) {
    x.resize(cols); y.resize(rows);
    for (Eigen::Index s = 0; s < cols; ++s)
        x[s] = integrateAxis(origin.sample + s, m.center.sample, m.sigmaPixels);
    for (Eigen::Index l = 0; l < rows; ++l)
        y[l] = integrateAxis(origin.line + l, m.center.line, m.sigmaPixels);
}

} // namespace fd::opnav::image::detail
