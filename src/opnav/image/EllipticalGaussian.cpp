#include "EllipticalGaussianInternal.hpp"

namespace fd::opnav::image {
EllipticalGaussianPixel evaluateEllipticalGaussianPixel(const EllipticalGaussian& m, PixelCoordinates pixel) {
    detail::validate(m);
    detail::validateOrigin(pixel);
    auto result = detail::integrateEllipse(m.center,m.heightDn,m.backgroundDn,detail::precision(m),pixel);
    const Eigen::Vector3d shape = result.partials.segment<3>(3);
    const double c = std::cos(m.angleRadians), s = std::sin(m.angleRadians);
    const double major = m.sigmaMajorPixels, minor = m.sigmaMinorPixels;
    result.partials[3] = shape.dot(Eigen::Vector3d{c*c,c*s,s*s})*(-2/(major*major*major));
    result.partials[4] = shape.dot(Eigen::Vector3d{s*s,-c*s,c*c})*(-2/(minor*minor*minor));
    result.partials[5] = shape.dot(Eigen::Vector3d{-2*c*s,c*c-s*s,2*c*s})*(1/(major*major)-1/(minor*minor));
    return result;
}

Image renderEllipticalGaussian(ImageSize size, const EllipticalGaussian& m, PixelCoordinates origin) {
    detail::validate(m);
    detail::validateOrigin(origin);
    if (size.width <= 0 || size.height <= 0) throw std::invalid_argument("Image dimensions must be positive.");
    const auto q = detail::precision(m);
    Image result(size.height,size.width);
    for (int y = 0; y < size.height; ++y)
        for (int x = 0; x < size.width; ++x)
            result(y,x) = detail::integrateEllipse(m.center,m.heightDn,m.backgroundDn,q,
                {origin.sample+x,origin.line+y}).dn;
    return result;
}
} // namespace fd::opnav::image
