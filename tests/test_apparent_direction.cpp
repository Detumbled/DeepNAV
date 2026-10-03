#include "dynamics/LinearStateProvider.hpp"
#include "opnav/CameraModel.hpp"
#include "opnav/core/ApparentDirection.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace fd::opnav;
using DynamicsState = fd::dynamics::CartesianState;
using fd::dynamics::LinearStateProvider;

namespace {
constexpr double c = 299792.458;

void near(double actual, double expected, double tolerance = 1e-10) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error("Numerical comparison failed.");
}

void vectorNear(const Eigen::Vector3d& actual, const Eigen::Vector3d& expected,
                double tolerance = 1e-9) {
    near((actual - expected).norm(), 0.0, tolerance);
}

template<class Exception = std::invalid_argument, class Function>
void rejects(Function function) {
    try { function(); }
    catch (const Exception&) { return; }
    throw std::runtime_error("Expected exception was not thrown.");
}

void lightTime() {
    const DynamicsState observer{};
    const LinearStateProvider stationary(TdbEpoch(0), {{0, 0, 3*c}, {0, 0, 0}}, "J2000", "SSB");
    const auto fixed = LightTimeSolver().solve(TdbEpoch(0), observer, stationary);
    near(fixed.lightTimeSeconds, 3);
    near(fixed.emissionEpoch.secondsPastJ2000, -3);
    near(fixed.residualSeconds, 0);
    near(fixed.iterations, 1);

    // Independent analytic solution: c*tau = initial range - v_target*tau.
    for (double velocity : {-12.0, 12.0}) {
        const LinearStateProvider moving(TdbEpoch(0), {{0, 0, 3*c}, {0, 0, velocity}}, "J2000", "SSB");
        const auto solution = LightTimeSolver().solve(TdbEpoch(0), observer, moving);
        const double tau = 3*c / (c + velocity);
        near(solution.lightTimeSeconds, tau);
        near(solution.geometricLineOfSightKm.z(), c*tau, 1e-5);
        vectorNear(solution.targetStateAtEmission.positionKm,
                   moving.stateAt(solution.emissionEpoch).positionKm);
        near(solution.receptionEpoch.secondsPastJ2000 - solution.emissionEpoch.secondsPastJ2000,
             solution.lightTimeSeconds);
        if (solution.residualSeconds > 1e-9) throw std::runtime_error("Excess light-time residual.");
        rejects<std::runtime_error>([&] {
            (void)LightTimeSolver({1e-15, 1}).solve(TdbEpoch(0), observer, moving);
        });
    }
}

void aberrationAndProjection() {
    const LinearStateProvider target(TdbEpoch(0), {{2*c, 0, 0}, {0, 0, 0}}, "J2000", "SSB");
    DynamicsState observer{};
    const auto uncorrected = core::computeApparentDirection(TdbEpoch(0), observer, target);
    vectorNear(uncorrected.vectorKm, {2*c, 0, 0});
    observer.velocityKmPerSec = {0, 30, 0};
    const auto apparent = core::computeApparentDirection(TdbEpoch(0), observer, target);
    vectorNear(apparent.vectorKm, {2*c, 60, 0});
    near(apparent.unitDirection.norm(), 1);
    near(apparent.unitDirection.y() / apparent.unitDirection.x(), 30/c);

    // Camera looks along inertial +X; camera +Y stays inertial +Y.
    Eigen::Matrix3d rotation;
    rotation << 0, 0, -1, 0, 1, 0, 1, 0, 0;
    const CameraAttitude attitude(rotation);
    vectorNear(attitude.toCamera({1, 2, 3}), {-3, 2, 1});
    vectorNear(CameraAttitude(Eigen::Matrix3d::Identity()).toCamera({1, 2, 3}), {1, 2, 3});
    vectorNear(attitude.matrix().transpose() * attitude.toCamera({1, 2, 3}), {1, 2, 3});

    const CameraModel camera({50, {1000, 500}, 20*Eigen::Matrix2d::Identity()});
    const auto pixel = camera.project(apparent.vectorKm, attitude);
    near(pixel.sample, 1000);
    near(pixel.line, 500 + 1000*30/c);
    const auto normalizedPixel = camera.project(apparent.unitDirection, attitude);
    near(normalizedPixel.sample, pixel.sample);
    near(normalizedPixel.line, pixel.line);
    rejects<std::domain_error>([&] { (void)camera.project(-apparent.vectorKm, attitude); });

    const Eigen::Vector3d offset{1e6, -2e6, 3e6};
    observer.positionKm += offset;
    const LinearStateProvider translated(TdbEpoch(0),
        {Eigen::Vector3d{2*c, 0, 0} + offset, {0, 0, 0}}, "J2000", "SSB");
    vectorNear(core::computeApparentDirection(TdbEpoch(0), observer, translated).vectorKm,
               apparent.vectorKm);
}

void movingGeometry() {
    const TdbEpoch reception(100);
    const DynamicsState observer{{100, -200, 300}, {3, -8, 4}};
    const Eigen::Vector3d relative{c, 0.5*c, 2*c};
    const Eigen::Vector3d velocity{10, 20, -5};
    const LinearStateProvider target(reception,
        {observer.positionKm + relative, velocity}, "J2000", "SSB");
    // Solve |relative - velocity*tau|^2 = c^2*tau^2 independently.
    const double a = c*c - velocity.squaredNorm();
    const double b = relative.dot(velocity);
    const double tau = (std::sqrt(b*b + a*relative.squaredNorm()) - b) / a;
    const auto result = core::computeApparentDirection(reception, observer, target, {1e-12, 8});
    near(result.lightTime.lightTimeSeconds, tau, 1e-12);
    vectorNear(result.lightTime.geometricLineOfSightKm, relative - velocity*tau, 1e-6);
    vectorNear(result.vectorKm, relative + (observer.velocityKmPerSec - velocity)*tau, 1e-6);
}

class InvalidProvider final : public fd::dynamics::StateProvider {
public:
    mutable std::size_t calls{0};
    DynamicsState stateAt(TdbEpoch) const override {
        ++calls;
        return {{0, 0, calls == 1 ? c : std::numeric_limits<double>::quiet_NaN()}, {0, 0, 0}};
    }
    std::string_view referenceFrame() const noexcept override { return "J2000"; }
    std::string_view referenceOrigin() const noexcept override { return "SSB"; }
};

void invalidInputs() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const LinearStateProvider target(TdbEpoch(0), {{0, 0, c}, {0, 0, 0}}, "J2000", "SSB");
    for (double invalid : {nan, inf}) {
        rejects([&] { (void)LightTimeSolver().solve(TdbEpoch(invalid), {}, target); });
        rejects([&] { (void)LightTimeSolver({invalid, 8}).solve(TdbEpoch(0), {}, target); });
        rejects([&] { (void)core::computeApparentDirection(TdbEpoch(0), {{}, {invalid, 0, 0}}, target); });
    }
    rejects([&] { (void)LightTimeSolver({0, 8}).solve(TdbEpoch(0), {}, target); });
    rejects([&] { (void)LightTimeSolver({1e-9, 0}).solve(TdbEpoch(0), {}, target); });
    rejects<std::runtime_error>([&] { (void)LightTimeSolver().solve(TdbEpoch(0), {}, InvalidProvider{}); });
    const LinearStateProvider coincident(TdbEpoch(0), {}, "J2000", "SSB");
    rejects<std::runtime_error>([&] { (void)LightTimeSolver().solve(TdbEpoch(0), {}, coincident); });

    Eigen::Matrix3d bad = Eigen::Matrix3d::Identity();
    bad(0, 0) = -1;
    rejects([&] { (void)CameraAttitude(bad); }); // Reflection is not an attitude.
    bad(0, 0) = 2;
    rejects([&] { (void)CameraAttitude(bad); });
    bad(0, 0) = nan;
    rejects([&] { (void)CameraAttitude(bad); });
    const CameraAttitude identity(Eigen::Matrix3d::Identity());
    rejects([&] { (void)identity.toCamera({0, 0, 0}); });
    rejects([&] { (void)identity.toCamera({nan, 0, 1}); });
}
} // namespace

int main() {
    try {
        lightTime();
        aberrationAndProjection();
        movingGeometry();
        invalidInputs();
        std::cout << "All apparent-direction, attitude, and projection pipeline tests passed.\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
