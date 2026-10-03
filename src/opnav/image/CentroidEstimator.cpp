#include "opnav/image/CentroidEstimator.hpp"
#include "CircularGaussianInternal.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <numbers>
#include <limits>

namespace fd::opnav::image {
namespace {
using Matrix5 = ParameterCovariance;

Parameters pack(const CircularGaussian& m) {
    return Parameters{m.center.sample, m.center.line, std::log(m.heightDn),
                      std::log(m.sigmaPixels), m.backgroundDn};
}
CircularGaussian unpack(const Parameters& p) {
    return {{p[0], p[1]}, std::exp(p[2]), std::exp(p[3]), p[4]};
}
bool inWindow(const CircularGaussian& m, const Image& dn, PixelCoordinates origin) {
    return m.center.sample >= origin.sample - 0.5
        && m.center.sample <= origin.sample + dn.cols() - 0.5
        && m.center.line >= origin.line - 0.5
        && m.center.line <= origin.line + dn.rows() - 0.5;
}

std::optional<CircularGaussian> initialize(const Image& dn, PixelCoordinates origin,
                                         const PixelMask* mask) {
    std::vector<double> border, all;
    for (Eigen::Index l = 0; l < dn.rows(); ++l)
        for (Eigen::Index s = 0; s < dn.cols(); ++s) {
            if (mask && (*mask)(l, s)) continue;
            all.push_back(dn(l, s));
            if (l == 0 || s == 0 || l == dn.rows()-1 || s == dn.cols()-1)
                border.push_back(dn(l, s));
        }
    auto& values = border.empty() ? all : border;
    auto middle = values.begin() + values.size()/2;
    std::nth_element(values.begin(), middle, values.end());
    const double background = *middle;
    double sum = 0, x = 0, y = 0;
    for (Eigen::Index l = 0; l < dn.rows(); ++l)
        for (Eigen::Index s = 0; s < dn.cols(); ++s) {
            if (mask && (*mask)(l, s)) continue;
            const double w = std::max(0.0, dn(l, s) - background);
            sum += w; x += w*s; y += w*l;
        }
    if (!(sum > 0) || !std::isfinite(sum)) return std::nullopt;
    x /= sum; y /= sum;
    double secondMoment = 0;
    for (Eigen::Index l = 0; l < dn.rows(); ++l)
        for (Eigen::Index s = 0; s < dn.cols(); ++s) {
            if (mask && (*mask)(l, s)) continue;
            const double w = std::max(0.0, dn(l, s) - background);
            secondMoment += w*((s-x)*(s-x) + (l-y)*(l-y));
        }
    const double sigma = std::clamp(std::sqrt(std::max(0.16, secondMoment/(2*sum) - 1.0/12)),
                                   0.4, 0.5*std::max(dn.rows(), dn.cols()));
    return CircularGaussian{{origin.sample+x, origin.line+y},
        sum/(2*std::numbers::pi*sigma*sigma), sigma, background};
}

struct NormalEquations {
    Matrix5 information{Matrix5::Zero()};
    Parameters gradient{Parameters::Zero()};
    double cost{0};
};

NormalEquations evaluate(const Image& dn, const Image& variance, const PixelMask* mask,
    PixelCoordinates origin, const CircularGaussian& m,
    std::vector<detail::AxisIntegral>& x, std::vector<detail::AxisIntegral>& y) {
    detail::fillAxes(m, dn.rows(), dn.cols(), origin, x, y);
    NormalEquations result;
    for (Eigen::Index l = 0; l < dn.rows(); ++l)
        for (Eigen::Index s = 0; s < dn.cols(); ++s) {
            if (mask && (*mask)(l, s)) continue;
            auto pixel = detail::combine(m, x[s], y[l]);
            // Positive height and sigma are optimized in logarithmic coordinates.
            pixel.partials[2] *= m.heightDn;
            pixel.partials[3] *= m.sigmaPixels;
            const double residual = dn(l, s) - pixel.dn;
            const double weight = 1 / variance(l, s);
            result.cost += residual*residual*weight;
            result.gradient.noalias() += (weight*residual)*pixel.partials;
            result.information.noalias() += weight*pixel.partials*pixel.partials.transpose();
        }
    if (!std::isfinite(result.cost) || !result.gradient.allFinite() || !result.information.allFinite())
        throw std::overflow_error("Gaussian normal equations overflowed.");
    return result;
}
} // namespace

GaussianFitResult fitCircularGaussian(const Image& dn, const Image& varianceDn2,
    PixelCoordinates origin, const PixelMask* mask, std::optional<CircularGaussian> initial,
    GaussianFitOptions options) {
    detail::validateOrigin(origin);
    if (dn.rows() < 3 || dn.cols() < 3 || dn.rows() != varianceDn2.rows()
        || dn.cols() != varianceDn2.cols()
        || (mask && (mask->rows() != dn.rows() || mask->cols() != dn.cols())))
        throw std::invalid_argument("Fit needs matching image/variance/mask grids, at least 3x3.");
    if (options.maxIterations == 0 || !std::isfinite(options.stepTolerance) || options.stepTolerance <= 0
        || !std::isfinite(options.minimumHeightSnr) || options.minimumHeightSnr < 0)
        throw std::invalid_argument("Invalid Gaussian fit options.");
    std::size_t usable = 0;
    for (Eigen::Index i = 0; i < dn.size(); ++i) {
        if (mask && mask->data()[i]) continue;
        if (!std::isfinite(dn.data()[i]) || !std::isfinite(varianceDn2.data()[i])
            || varianceDn2.data()[i] <= 0 || !std::isfinite(1/varianceDn2.data()[i]))
            throw std::invalid_argument("Unmasked pixels need finite DN and positive finite variance.");
        ++usable;
    }
    if (usable <= 5) throw std::invalid_argument("Five-parameter fit needs more than five usable pixels.");
    GaussianFitResult result;
    result.degreesOfFreedom = usable - 5;
    if (!initial) initial = initialize(dn, origin, mask);
    if (!initial) { result.status = FitStatus::NoSignal; return result; }
    detail::validate(*initial);
    if (!inWindow(*initial, dn, origin))
        throw std::invalid_argument("Initial center must lie inside the fitting window.");
    Parameters p = pack(*initial);
    std::vector<detail::AxisIntegral> x, y; // Reuse axis storage across iterations.
    auto normal = evaluate(dn, varianceDn2, mask, origin, *initial, x, y);
    double damping = 1e-3;
    for (std::size_t iteration = 0; iteration < options.maxIterations; ++iteration) {
        result.iterations = iteration + 1;
        if ((normal.information.diagonal().array() <= 0).any()) {
            result.status = FitStatus::Singular; break;
        }
        const Parameters scale = normal.information.diagonal().array().sqrt().inverse();
        const Matrix5 scaled = scale.asDiagonal()*normal.information*scale.asDiagonal();
        const Parameters gradient = scale.asDiagonal()*normal.gradient;
        if (gradient.norm() <= options.stepTolerance) { result.status = FitStatus::Converged; break; }
        bool accepted = false;
        for (int trial = 0; trial < 16; ++trial) {
            Matrix5 damped = scaled;
            damped.diagonal().array() += damping;
            const Parameters step = scale.asDiagonal()*damped.ldlt().solve(gradient);
            const Parameters candidateP = p + step;
            const auto candidate = unpack(candidateP);
            if (candidateP.allFinite() && std::isfinite(candidate.heightDn) && candidate.heightDn > 0
                && std::isfinite(candidate.sigmaPixels) && candidate.sigmaPixels >= 0.05
                && candidate.sigmaPixels <= 2*std::max(dn.rows(), dn.cols())
                && inWindow(candidate, dn, origin)) {
                const auto next = evaluate(dn, varianceDn2, mask, origin, candidate, x, y);
                if (next.cost <= normal.cost) {
                    const double relativeBackgroundStep = std::abs(step[4])/(1+std::abs(p[4]));
                    p = candidateP; normal = next; accepted = true;
                    damping = std::max(1e-12, damping/3);
                    if (step.head<4>().cwiseAbs().maxCoeff() <= options.stepTolerance
                        && relativeBackgroundStep <= options.stepTolerance)
                        result.status = FitStatus::Converged;
                    break;
                }
            }
            damping *= 10;
        }
        if (!accepted || result.status == FitStatus::Converged) break;
    }
    result.model = unpack(p);
    result.chiSquared = normal.cost;
    if (result.status != FitStatus::Converged) return result;

    const Parameters scale = normal.information.diagonal().array().sqrt().inverse();
    const Matrix5 scaled = scale.asDiagonal()*normal.information*scale.asDiagonal();
    const Eigen::SelfAdjointEigenSolver<Matrix5> eigenvalues(scaled, Eigen::EigenvaluesOnly);
    if (eigenvalues.info() != Eigen::Success || !scale.allFinite()
        || eigenvalues.eigenvalues()[0] <= 1e-10*eigenvalues.eigenvalues()[4]) {
        result.status = FitStatus::Singular; return result;
    }
    Parameters physicalScale = scale;
    physicalScale[2] *= result.model.heightDn;
    physicalScale[3] *= result.model.sigmaPixels;
    Matrix5 covariance = physicalScale.asDiagonal()*scaled.ldlt().solve(Matrix5::Identity())
        *physicalScale.asDiagonal();
    covariance = (0.5*covariance + 0.5*covariance.transpose()).eval();
    if (!covariance.allFinite() || (covariance.diagonal().array() <= 0).any()) {
        result.status = FitStatus::Singular; return result;
    }
    if (result.model.heightDn / std::sqrt(covariance(2, 2)) < options.minimumHeightSnr) {
        result.status = FitStatus::NoSignal; return result;
    }
    result.parameterCovariance = covariance;
    result.measurement = PixelMeasurement{result.model.center, covariance.topLeftCorner<2, 2>()};
    return result;
}

GaussianFitResult fitBrightestCircularGaussian(const Image& dn, const Image& varianceDn2,
    int windowRadius, PixelCoordinates origin, const PixelMask* mask, GaussianFitOptions options) {
    detail::validateOrigin(origin);
    if (windowRadius < 1 || dn.rows() < 3 || dn.cols() < 3
        || dn.rows() != varianceDn2.rows() || dn.cols() != varianceDn2.cols()
        || (mask && (mask->rows() != dn.rows() || mask->cols() != dn.cols())))
        throw std::invalid_argument("Source search needs matching grids, at least 3x3, and positive radius.");
    Eigen::Index peakSample = 0, peakLine = 0;
    double peak = -std::numeric_limits<double>::infinity();
    for (Eigen::Index l = 0; l < dn.rows(); ++l)
        for (Eigen::Index s = 0; s < dn.cols(); ++s) {
            if (mask && (*mask)(l, s)) continue;
            if (!std::isfinite(dn(l, s)) || !std::isfinite(varianceDn2(l, s))
                || varianceDn2(l, s) <= 0 || !std::isfinite(1/varianceDn2(l, s)))
                throw std::invalid_argument("Unmasked pixels need finite DN and positive finite variance.");
            if (dn(l, s) > peak) { peak = dn(l, s); peakSample = s; peakLine = l; }
        }
    if (!std::isfinite(peak)) throw std::invalid_argument("Source search has no usable pixels.");
    const Eigen::Index diameter = 2*static_cast<Eigen::Index>(windowRadius) + 1;
    const Eigen::Index width = std::min(diameter, dn.cols());
    const Eigen::Index height = std::min(diameter, dn.rows());
    // Shift edge windows inward to retain as much background as possible.
    const Eigen::Index firstSample = std::clamp(peakSample-windowRadius, Eigen::Index{0}, dn.cols()-width);
    const Eigen::Index firstLine = std::clamp(peakLine-windowRadius, Eigen::Index{0}, dn.rows()-height);
    const Image window = dn.block(firstLine, firstSample, height, width);
    const Image variance = varianceDn2.block(firstLine, firstSample, height, width);
    PixelMask windowMask;
    if (mask) windowMask = mask->block(firstLine, firstSample, height, width);
    return fitCircularGaussian(window, variance,
        {origin.sample+firstSample, origin.line+firstLine}, mask ? &windowMask : nullptr,
        std::nullopt, options);
}

} // namespace fd::opnav::image
