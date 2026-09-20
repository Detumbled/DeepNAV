#include "opnav/CameraModel.hpp"
#include "opnav/DistortionModel.hpp"

#include <Eigen/Core>

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

bool approximatelyEqual(
    double actual,
    double expected,
    double tolerance = 1.0e-12) {

    return std::abs(actual - expected) <= tolerance;
}

void requirePixel(
    const fd::opnav::PixelCoordinates& actual,
    double expectedSample,
    double expectedLine,
    double tolerance = 1.0e-12) {

    if (!approximatelyEqual(
            actual.sample, expectedSample, tolerance) ||
        !approximatelyEqual(
            actual.line, expectedLine, tolerance)) {

        std::cerr
            << "Expected pixel: ("
            << expectedSample << ", "
            << expectedLine << ")\n"
            << "Actual pixel:   ("
            << actual.sample << ", "
            << actual.line << ")\n";

        std::abort();
    }
}

fd::opnav::CameraIntrinsics makeTestIntrinsics() {
    Eigen::Matrix2d focalPlaneToPixel;

    focalPlaneToPixel <<
        20.0,  0.0,
         0.0, 20.0;

    return fd::opnav::CameraIntrinsics {
        .focalLength = 50.0,
        .principalPoint = {
            .sample = 1000.0,
            .line   = 500.0
        },
        .focalPlaneToPixelMatrix_PixelPermm =
            focalPlaneToPixel
    };
}

void testNoDistortion() {
    using namespace fd::opnav;

    const CameraModel camera {makeTestIntrinsics()};

    /*
     * Normalized coordinates:
     *
     * xn = 0.1 / 1.0 = 0.1
     * yn = 0.2 / 1.0 = 0.2
     *
     * Focal-plane coordinates:
     *
     * x = 50 * 0.1 = 5 mm
     * y = 50 * 0.2 = 10 mm
     *
     * Pixel offsets:
     *
     * ds = 20 * 5  = 100 px
     * dl = 20 * 10 = 200 px
     *
     * Final pixel:
     *
     * s = 1000 + 100 = 1100
     * l =  500 + 200 = 700
     */

    const PixelCoordinates pixel =
        camera.project(Eigen::Vector3d {0.1, 0.2, 1.0});

    requirePixel(pixel, 1100.0, 700.0);

    std::cout << "PASS: no distortion\n";
}

void testZeroOpenCvDistortion() {
    using namespace fd::opnav;

    const auto distortion =
        std::make_shared<OpenCVDistortionModel>(
            OpenCVDistortionCoefficients {});

    const CameraModel camera {
        makeTestIntrinsics(),
        distortion
    };

    const PixelCoordinates pixel =
        camera.project(Eigen::Vector3d {0.1, 0.2, 1.0});

    requirePixel(pixel, 1100.0, 700.0);

    std::cout << "PASS: zero OpenCV coefficients\n";
}

void testRadialK1Distortion() {
    using namespace fd::opnav;

    OpenCVDistortionCoefficients coefficients;
    coefficients.K1 = 0.1;

    const auto distortion =
        std::make_shared<OpenCVDistortionModel>(
            coefficients);

    const CameraModel camera {
        makeTestIntrinsics(),
        distortion
    };

    /*
     * xn = 0.1
     * yn = 0.2
     *
     * r² = 0.1² + 0.2² = 0.05
     *
     * radial scale:
     *
     * 1 + K1*r² = 1 + 0.1*0.05 = 1.005
     *
     * distorted normalized coordinates:
     *
     * xd = 0.1 * 1.005 = 0.1005
     * yd = 0.2 * 1.005 = 0.201
     *
     * focal plane:
     *
     * xf = 50 * 0.1005 = 5.025 mm
     * yf = 50 * 0.201  = 10.05 mm
     *
     * pixels:
     *
     * s = 1000 + 20*5.025 = 1100.5
     * l =  500 + 20*10.05 = 701.0
     */

    const PixelCoordinates pixel =
        camera.project(Eigen::Vector3d {0.1, 0.2, 1.0});

    requirePixel(pixel, 1100.5, 701.0, 1.0e-10);

    std::cout << "PASS: radial K1 distortion\n";
}

void testBoresight() {
    using namespace fd::opnav;

    OpenCVDistortionCoefficients coefficients;
    coefficients.K1 = 0.5;
    coefficients.P1 = 0.1;
    coefficients.P2 = -0.1;

    const CameraModel camera {
        makeTestIntrinsics(),
        std::make_shared<OpenCVDistortionModel>(
            coefficients)
    };

    // At r = 0, all distortion contributions vanish.
    const PixelCoordinates pixel =
        camera.project(Eigen::Vector3d {0.0, 0.0, 1.0});

    requirePixel(pixel, 1000.0, 500.0);

    std::cout << "PASS: distorted boresight\n";
}

void testTargetBehindCamera() {
    using namespace fd::opnav;

    const CameraModel camera {makeTestIntrinsics()};

    bool exceptionThrown = false;

    try {
        [[maybe_unused]]
        const auto pixel =
            camera.project(Eigen::Vector3d {0.1, 0.2, -1.0});
    } catch (const std::domain_error&) {
        exceptionThrown = true;
    }

    assert(exceptionThrown);

    std::cout << "PASS: target behind camera rejected\n";
}

} // namespace

int main() {
    testNoDistortion();
    testZeroOpenCvDistortion();
    testRadialK1Distortion();
    testBoresight();
    testTargetBehindCamera();

    std::cout << "All optical camera tests passed.\n";

    return 0;
}