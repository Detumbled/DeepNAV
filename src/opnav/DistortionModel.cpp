#include "opnav/DistortionModel.hpp"
#include "opnav/Types.hpp"

#include <cmath>
#include <stdexcept>

namespace fd::opnav {
            FocalPlaneCoordinates OpenCVDistortionModel::distort(
            const FocalPlaneCoordinates& ideal) const {

            const double x = ideal.x();
            const double y = ideal.y();

            const double r2 = x * x + y * y;
            const double r4 = r2 * r2;
            const double r6 = r4 * r2;

            const double numerator =
                1.0
                + coefficients_.K1 * r2
                + coefficients_.K2 * r4
                + coefficients_.K3 * r6;

            const double denominator =
                1.0
                + coefficients_.K4 * r2
                + coefficients_.K5 * r4
                + coefficients_.K6 * r6;

            if (!std::isfinite(denominator) ||
                std::abs(denominator) < 1.0e-14) {
                throw std::domain_error(
                    "OpenCV distortion denominator is singular.");
            }

            const double radialScale = numerator / denominator;

            const double tangentialX =
                2.0 * coefficients_.P1 * x * y
                + coefficients_.P2 * (r2 + 2.0 * x * x);

            const double tangentialY =
                coefficients_.P1 * (r2 + 2.0 * y * y)
                + 2.0 * coefficients_.P2 * x * y;

            return FocalPlaneCoordinates{
                 x * radialScale + tangentialX,
                 y * radialScale + tangentialY
            };
        }



}