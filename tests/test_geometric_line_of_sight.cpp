#include "opnav/core/GeometricLineOfSight.hpp"

#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

constexpr double kTolerance = 1.0e-12;

void requireNear(
    double actual,
    double expected,
    double tolerance = kTolerance) {

    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(
            "Expected " + std::to_string(expected)
            + ", obtained " + std::to_string(actual));
    }
}

void requireVectorNear(
    const Eigen::Vector3d& actual,
    const Eigen::Vector3d& expected,
    double tolerance = kTolerance) {

    if ((actual - expected).norm() > tolerance) {
        throw std::runtime_error("Vector comparison failed.");
    }
}

void testKnownGeometry() {
    const fd::dynamics::CartesianState camera {
        .positionKm = {1000.0, -2000.0, 500.0},
        .velocityKmPerSec = {1.0, 2.0, 3.0}
    };

    const fd::dynamics::CartesianState target {
        .positionKm = {4000.0, 2000.0, 12500.0},
        .velocityKmPerSec = {-4.0, 5.0, 6.0}
    };

    const auto result =
        fd::opnav::core::computeGeometricLineOfSight(camera, target);

    requireVectorNear(
        result.vectorKm,
        Eigen::Vector3d {3000.0, 4000.0, 12000.0});

    requireNear(result.rangeKm, 13000.0);

    requireVectorNear(
        result.unitDirection,
        Eigen::Vector3d {3.0 / 13.0, 4.0 / 13.0, 12.0 / 13.0});

    requireNear(result.unitDirection.norm(), 1.0);

    std::cout << "PASS: known geometric line of sight\n";
}

void testTranslationInvariance() {
    const Eigen::Vector3d translation {1.0e6, -2.0e6, 3.0e6};

    const fd::dynamics::CartesianState camera {
        .positionKm = {1000.0, -2000.0, 500.0},
        .velocityKmPerSec = {0.0, 0.0, 0.0}
    };

    const fd::dynamics::CartesianState target {
        .positionKm = {4000.0, 2000.0, 12500.0},
        .velocityKmPerSec = {0.0, 0.0, 0.0}
    };

    auto translatedCamera = camera;
    auto translatedTarget = target;
    translatedCamera.positionKm += translation;
    translatedTarget.positionKm += translation;

    const auto original =
        fd::opnav::core::computeGeometricLineOfSight(camera, target);

    const auto translated =
        fd::opnav::core::computeGeometricLineOfSight(
            translatedCamera,
            translatedTarget);

    requireVectorNear(translated.vectorKm, original.vectorKm);
    requireVectorNear(translated.unitDirection, original.unitDirection);
    requireNear(translated.rangeKm, original.rangeKm);

    std::cout << "PASS: translation invariance\n";
}

void testCoincidentPositionsRejected() {
    const fd::dynamics::CartesianState camera {
        .positionKm = {100.0, 200.0, 300.0},
        .velocityKmPerSec = {0.0, 0.0, 0.0}
    };

    const fd::dynamics::CartesianState target {
        .positionKm = {100.0, 200.0, 300.0},
        .velocityKmPerSec = {1.0, 0.0, 0.0}
    };

    bool exceptionThrown = false;

    try {
        [[maybe_unused]] const auto result =
            fd::opnav::core::computeGeometricLineOfSight(camera, target);
    } catch (const std::domain_error&) {
        exceptionThrown = true;
    }

    if (!exceptionThrown) {
        throw std::runtime_error(
            "Coincident positions were not rejected.");
    }

    std::cout << "PASS: coincident positions rejected\n";
}

void testNonFiniteStateRejected() {
    fd::dynamics::CartesianState camera {
        .positionKm = {0.0, 0.0, 0.0},
        .velocityKmPerSec = {0.0, 0.0, 0.0}
    };

    const fd::dynamics::CartesianState target {
        .positionKm = {1000.0, 0.0, 0.0},
        .velocityKmPerSec = {0.0, 0.0, 0.0}
    };

    camera.positionKm.x() =
        std::numeric_limits<double>::quiet_NaN();

    bool exceptionThrown = false;

    try {
        [[maybe_unused]] const auto result =
            fd::opnav::core::computeGeometricLineOfSight(camera, target);
    } catch (const std::invalid_argument&) {
        exceptionThrown = true;
    }

    if (!exceptionThrown) {
        throw std::runtime_error(
            "Non-finite state was not rejected.");
    }

    std::cout << "PASS: non-finite state rejected\n";
}

} // namespace

int main() {
    try {
        testKnownGeometry();
        testTranslationInvariance();
        testCoincidentPositionsRejected();
        testNonFiniteStateRejected();

        std::cout
            << "All geometric line-of-sight tests passed.\n";

        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "FAIL: " << exception.what() << '\n';
        return 1;
    }
}
