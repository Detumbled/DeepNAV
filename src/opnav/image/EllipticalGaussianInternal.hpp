#pragma once

#include "opnav/image/EllipticalGaussian.hpp"
#include "CircularGaussianInternal.hpp"
#include <array>

namespace fd::opnav::image::detail {

inline void validate(const EllipticalGaussian& m) {
    validateOrigin(m.center);
    if (!std::isfinite(m.heightDn) || m.heightDn <= 0 ||
        !std::isfinite(m.sigmaMajorPixels) || m.sigmaMajorPixels < 0.2 ||
        !std::isfinite(m.sigmaMinorPixels) || m.sigmaMinorPixels < 0.2 ||
        !std::isfinite(m.angleRadians) || !std::isfinite(m.backgroundDn))
        throw std::invalid_argument("Elliptical Gaussian needs finite parameters, positive height, and widths >= 0.2 px.");
}

inline Eigen::Matrix2d precision(const EllipticalGaussian& m) {
    const double c = std::cos(m.angleRadians), s = std::sin(m.angleRadians);
    Eigen::Matrix2d rotation;
    rotation << c,-s,s,c;
    return rotation * Eigen::Vector2d{1/(m.sigmaMajorPixels*m.sigmaMajorPixels),
        1/(m.sigmaMinorPixels*m.sigmaMinorPixels)}.asDiagonal() * rotation.transpose();
}

struct Quadrature { std::array<double,16> nodes, weights; };
inline const Quadrature& quadrature() {
    static const Quadrature rule = [] {
        Quadrature q{};
        for (int i = 0; i < 8; ++i) {
            double z = std::cos(std::numbers::pi*(i+0.75)/16.5), derivative = 0;
            for (int iteration = 0; iteration < 30; ++iteration) {
                double p = 1, previous = 0;
                for (int n = 1; n <= 16; ++n) {
                    const double next = ((2*n-1)*z*p-(n-1)*previous)/n;
                    previous = p; p = next;
                }
                derivative = 16*(z*p-previous)/(z*z-1);
                const double step = p/derivative;
                z -= step;
                if (std::abs(step)<1e-15) break;
            }
            // Map [-1,1] to a pixel's [-0.5,0.5] offsets and weights.
            q.nodes[i] = -z/2; q.nodes[15-i] = z/2;
            q.weights[i] = q.weights[15-i] = 1/((1-z*z)*derivative*derivative);
        }
        return q;
    }();
    return rule;
}

// Derivatives in [sample,line,height,Qxx,Qxy,Qyy,background], Q = precision.
// All model-dependent coefficients are computed once per rendering/fit iteration.
inline EllipticalGaussianPixel integrateEllipse(PixelCoordinates center, double height,
    double background, const Eigen::Matrix2d& q, PixelCoordinates pixel) {
    const auto& rule = quadrature();
    EllipticalGaussianPixel result{background,EllipticalParameters::Zero()};
    result.partials[6] = 1;
    const double x = pixel.sample-center.sample, y = pixel.line-center.line;
    for (int j = 0; j < 16; ++j) {
        const double dy = y+rule.nodes[j];
        for (int i = 0; i < 16; ++i) {
            const double dx = x+rule.nodes[i];
            const double qx = q(0,0)*dx+q(0,1)*dy, qy = q(0,1)*dx+q(1,1)*dy;
            const double shape = rule.weights[i]*rule.weights[j]*std::exp(-0.5*(dx*qx+dy*qy));
            const double signal = height*shape;
            result.dn += signal;
            result.partials[0] += signal*qx;
            result.partials[1] += signal*qy;
            result.partials[2] += shape;
            result.partials[3] -= 0.5*signal*dx*dx;
            result.partials[4] -= signal*dx*dy;
            result.partials[5] -= 0.5*signal*dy*dy;
        }
    }
    if (!std::isfinite(result.dn) || !result.partials.allFinite())
        throw std::overflow_error("Elliptical Gaussian pixel calculation overflowed.");
    return result;
}
} // namespace fd::opnav::image::detail
