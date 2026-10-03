#include "opnav/image/CircularGaussian.hpp"
#include "CircularGaussianInternal.hpp"
#include <random>

namespace fd::opnav::image {

GaussianPixel evaluateGaussianPixel(const CircularGaussian& model, PixelCoordinates pixel) {
    detail::validate(model);
    detail::validateOrigin(pixel);
    return detail::combine(model,
        detail::integrateAxis(pixel.sample, model.center.sample, model.sigmaPixels),
        detail::integrateAxis(pixel.line, model.center.line, model.sigmaPixels));
}

Image renderCircularGaussian(ImageSize size, const CircularGaussian& model, PixelCoordinates origin) {
    detail::validate(model);
    detail::validateOrigin(origin);
    if (size.width <= 0 || size.height <= 0)
        throw std::invalid_argument("Image dimensions must be positive.");
    Image result(size.height, size.width);
    std::vector<detail::AxisIntegral> x, y;
    detail::fillAxes(model, size.height, size.width, origin, x, y);
    for (int l = 0; l < size.height; ++l)
        for (int s = 0; s < size.width; ++s)
            result(l, s) = detail::combine(model, x[s], y[l]).dn;
    return result;
}

SyntheticImage simulateCircularGaussian(ImageSize size, const CircularGaussian& model,
    DetectorNoise noise, std::uint64_t seed, PixelCoordinates origin) {
    if (!std::isfinite(noise.electronsPerDn) || noise.electronsPerDn <= 0
        || !std::isfinite(noise.readNoiseDn) || noise.readNoiseDn < 0
        || model.backgroundDn < 0)
        throw std::invalid_argument("Simulation needs positive gain and nonnegative background/read noise.");
    SyntheticImage result{renderCircularGaussian(size, model, origin), Image(size.height, size.width)};
    std::mt19937_64 generator(seed);
    std::normal_distribution<double> readNoise(0, 1);
    for (Eigen::Index i = 0; i < result.dn.size(); ++i) {
        const double mean = result.dn.data()[i];
        const double electrons = mean * noise.electronsPerDn;
        const double variance = mean / noise.electronsPerDn + noise.readNoiseDn*noise.readNoiseDn;
        if (!std::isfinite(electrons) || electrons > 1e15 || !std::isfinite(variance))
            throw std::overflow_error("Detector simulation exceeds supported count range.");
        result.varianceDn2.data()[i] = variance;
        const double shot = electrons > 0
            ? static_cast<double>(std::poisson_distribution<std::uint64_t>(electrons)(generator)) : 0;
        result.dn.data()[i] = shot / noise.electronsPerDn + noise.readNoiseDn*readNoise(generator);
    }
    return result;
}

} // namespace fd::opnav::image
