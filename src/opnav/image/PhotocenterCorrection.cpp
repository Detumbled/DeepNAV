#include "opnav/image/PhotocenterCorrection.hpp"
#include "CircularGaussianInternal.hpp"

namespace fd::opnav::image {
namespace {
void validateCovariance(const Eigen::Matrix2d& covariance) {
    if (!covariance.allFinite() || !covariance.isApprox(covariance.transpose(), 1e-12)
        || covariance(0, 0) < 0 || covariance(1, 1) < 0
        || std::abs(covariance(0, 1)) > std::sqrt(covariance(0, 0))*std::sqrt(covariance(1, 1)))
        throw std::invalid_argument("Pixel covariance must be finite, symmetric and positive semidefinite.");
}
}

PhotocenterOffset brightnessPhotocenterOffset(const Image& brightness,
    PixelCoordinates geometricCenter, PixelCoordinates origin) {
    detail::validateOrigin(origin);
    detail::validateOrigin(geometricCenter);
    if (brightness.size() == 0 || !brightness.allFinite() || (brightness.array() < 0).any())
        throw std::invalid_argument("Brightness model must be finite, nonnegative and nonempty.");
    // Scale before summing so arbitrary brightness normalization cancels safely.
    const double peak = brightness.maxCoeff();
    if (!(peak > 0)) throw std::invalid_argument("Brightness model has no illuminated pixels.");
    double sum = 0;
    Eigen::Vector2d moment = Eigen::Vector2d::Zero();
    for (Eigen::Index l = 0; l < brightness.rows(); ++l)
        for (Eigen::Index s = 0; s < brightness.cols(); ++s) {
            const double weight = brightness(l, s)/peak;
            sum += weight;
            moment += weight*Eigen::Vector2d(origin.sample+s-geometricCenter.sample,
                                            origin.line+l-geometricCenter.line);
        }
    PhotocenterOffset result;
    result.pixels = moment/sum;
    if (!result.pixels.allFinite()) throw std::overflow_error("Photocenter moment overflowed.");
    return result;
}

PixelMeasurement correctPhotocenter(const PixelMeasurement& photocenter, const PhotocenterOffset& offset) {
    detail::validateOrigin(photocenter.center);
    if (!offset.pixels.allFinite()) throw std::invalid_argument("Photocenter offset must be finite.");
    validateCovariance(photocenter.covariance);
    validateCovariance(offset.covariance);
    PixelMeasurement result{{photocenter.center.sample-offset.pixels.x(),
                             photocenter.center.line-offset.pixels.y()},
                            photocenter.covariance+offset.covariance};
    detail::validateOrigin(result.center);
    if (!result.covariance.allFinite()) throw std::overflow_error("Corrected covariance overflowed.");
    return result;
}

} // namespace fd::opnav::image
