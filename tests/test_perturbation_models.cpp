#include "filters/CartesianPropagator.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/SRP.hpp"
#include "perturbations/Eclipse.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace fd::perturbations;
namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class DerivedA, class DerivedB>
void relativeNear(const Eigen::MatrixBase<DerivedA>& actual,
                  const Eigen::MatrixBase<DerivedB>& expected, double tolerance) {
    if (expected.norm() > 0 && (actual - expected).norm() / expected.norm() >= tolerance)
        std::cerr << "Relative comparison error: " << (actual - expected).norm() / expected.norm()
                  << "; tolerance: " << tolerance << '\n';
    require(actual.allFinite() && expected.norm() > 0.0 &&
                (actual - expected).norm() / expected.norm() < tolerance,
            "Force/partial comparison failed.");
}
template <class Action> void rejects(Action action) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error("Expected invalid force input to be rejected.");
}
void checkPartials(const AccelerationFunction& model, const Eigen::Vector3d& position,
                   double step) {
    Eigen::Matrix3d numeric;
    for (int j = 0; j < 3; ++j) {
        Eigen::Vector3d plus = position, minus = position;
        plus[j] += step;
        minus[j] -= step;
        numeric.col(j) = (model(0, plus, Eigen::Vector3d::Zero()).acceleration -
                          model(0, minus, Eigen::Vector3d::Zero()).acceleration) /
                         (2 * step);
    }
    relativeNear(model(0, position, Eigen::Vector3d::Zero()).positionJacobian, numeric, 1e-6);
}

void forceLawsAndPartials() {
    const SolarRadiationPressure srp("EARTH", "J2000", 1.3, 20, 1000);
    const Eigen::Vector3d sun(SolarRadiationPressure::kAu_km, 0, 0), zero = Eigen::Vector3d::Zero();
    const double expected = SolarRadiationPressure::kSolarPressure1AU_N_m2 * 1.3 * 20 / 1000 / 1000;
    relativeNear(srp.evaluateAtSunPosition(zero, sun).acceleration,
                 Eigen::Vector3d(-expected, 0, 0), 1e-14);
    relativeNear(srp.evaluateAtSunPosition(-sun, sun).acceleration,
                 Eigen::Vector3d(-expected / 4, 0, 0), 1e-14);
    const Eigen::Vector3d position(20000, 4000, 3000);
    const auto radiation = [srp, sun](double, const Eigen::Vector3d& r, const Eigen::Vector3d&) {
        return srp.evaluateAtSunPosition(r, sun);
    };
    checkPartials(radiation, position, 100);
    const auto gravity = pointMassGravity(398600.4418);
    checkPartials(gravity, position, .1);
    const Eigen::Vector3d moon(300000, 200000, -100000);
    const auto tidal = [moon](double, const Eigen::Vector3d& r, const Eigen::Vector3d&) {
        return thirdBodyGravity(4902.800066, r, moon);
    };
    const Eigen::Vector3d distance = moon - position;
    const Eigen::Vector3d reference =
        4902.800066 * (distance / std::pow(distance.norm(), 3) - moon / std::pow(moon.norm(), 3));
    relativeNear(tidal(0, position, zero).acceleration, reference, 1e-12);
    checkPartials(tidal, position, 1);
    const auto combined = sumAccelerations({gravity, tidal, radiation})(0, position, zero);
    relativeNear(combined.acceleration,
                 gravity(0, position, zero).acceleration + tidal(0, position, zero).acceleration +
                     radiation(0, position, zero).acceleration,
                 1e-14);
    relativeNear(combined.positionJacobian,
                 gravity(0, position, zero).positionJacobian +
                     tidal(0, position, zero).positionJacobian +
                     radiation(0, position, zero).positionJacobian,
                 1e-14);
    rejects([&] { (void)srp.evaluateAtSunPosition(sun, sun); });
    rejects([&] { (void)pointMassGravity(-1); });
    rejects([&] { (void)thirdBodyGravity(4902, position, zero); });
    rejects([&] { (void)sumAccelerations({AccelerationFunction{}}); });
}

void propagatedPartialsAndNoise() {
    const Eigen::Vector3d sun(SolarRadiationPressure::kAu_km, 0, 0);
    const SolarRadiationPressure srp("EARTH", "J2000", 1.3, 20, 1000);
    const auto gravity = pointMassGravity(398600.4418);
    const auto forces = sumAccelerations(
        {gravity, [sun, srp](double, const Eigen::Vector3d& r, const Eigen::Vector3d&) {
             return srp.evaluateAtSunPosition(r, sun);
         }});
    fd::filters::CartesianPropagationConfig config;
    config.integrator.absoluteTolerance = 1e-12;
    config.integrator.relativeTolerance = 1e-12;
    config.integrator.maximumStep = 30;
    config.accelerationDiffusion = 1e-12 * Eigen::Vector3d(1, 2, 3).asDiagonal();
    const fd::filters::CartesianPropagator propagate(forces, config);
    Eigen::VectorXd initial(6);
    initial << 20000, 4000, 3000, -1, 4, 1;
    const auto result = propagate(0, 600, initial);
    Eigen::Matrix<double, 6, 6> numeric;
    for (int j = 0; j < 6; ++j) {
        const double step = j < 3 ? .1 : 1e-4;
        Eigen::VectorXd plus = initial, minus = initial;
        plus[j] += step;
        minus[j] -= step;
        numeric.col(j) =
            (propagate(0, 600, plus).state - propagate(0, 600, minus).state) / (2 * step);
    }
    relativeNear(result.transition, numeric, 1e-7);
    const auto first = propagate(0, 240, initial), second = propagate(240, 600, first.state);
    relativeNear(result.processCovariance,
                 second.transition * first.processCovariance * second.transition.transpose() +
                     second.processCovariance,
                 1e-8);
    relativeNear(result.state, second.state, 1e-11);
    const auto unperturbed = fd::filters::CartesianPropagator(gravity, config)(0, 600, initial);
    require((result.state - unperturbed.state).head<3>().norm() > 1e-6,
            "SRP had no propagated effect.");
}

void eclipseGeometryAndPropagation() {
    const Eigen::Vector3d sun(SolarRadiationPressure::kAu_km, 0, 0);
    constexpr double earthRadius = 6378.137, sunRadius = 695700;
    const auto lit = evaluateEclipse({20000, 0, 0}, sun, earthRadius, sunRadius);
    const auto dark = evaluateEclipse({-20000, 0, 0}, sun, earthRadius, sunRadius);
    require(lit.illumination == 1 && lit.flag == EclipseFlag::Sunlit &&
                dark.illumination == 0 && dark.flag == EclipseFlag::Umbra,
            "Sunlit/umbra classification failed.");
    double low = 6200, high = 6550;
    for (int i = 0; i < 40; ++i) {
        const double y = (low + high) / 2;
        if (evaluateEclipse({-20000, y, 0}, sun, earthRadius, sunRadius).illumination < .5)
            low = y;
        else
            high = y;
    }
    const Eigen::Vector3d edge(-20000, (low + high) / 2, 0);
    const auto partial = evaluateEclipse(edge, sun, earthRadius, sunRadius);
    require(partial.flag == EclipseFlag::Penumbra && std::abs(partial.illumination - .5) < 1e-8,
            "Penumbra fraction failed.");
    // Independent ray/sphere intersection across a uniformly sampled angular solar disk.
    const Eigen::Vector3d axis = (sun - edge).normalized();
    const Eigen::Vector3d u = axis.unitOrthogonal(), v = axis.cross(u);
    const double angularRadius = std::asin(sunRadius / (sun - edge).norm());
    int illuminated = 0, total = 0;
    for (int i = -80; i <= 80; ++i)
        for (int j = -80; j <= 80; ++j) {
            const double x = i / 80.0, y = j / 80.0, diskRadius = std::hypot(x, y);
            if (diskRadius > 1)
                continue;
            Eigen::Vector3d ray = axis;
            if (diskRadius > 0)
                ray = std::cos(angularRadius * diskRadius) * axis +
                      std::sin(angularRadius * diskRadius) * (x * u + y * v) / diskRadius;
            const double distanceAlongRay = -edge.dot(ray);
            const bool blocked = distanceAlongRay > 0 &&
                                 edge.squaredNorm() - distanceAlongRay * distanceAlongRay <
                                     earthRadius * earthRadius;
            illuminated += !blocked;
            ++total;
        }
    require(std::abs(partial.illumination - double(illuminated) / total) < .01,
            "Angular disk overlap disagrees with ray/sphere shadow geometry.");
    const SolarRadiationPressure srp("EARTH", "J2000", 1.3, 20, 1000);
    const auto radiation = [=](double, const Eigen::Vector3d& r, const Eigen::Vector3d&) {
        return srp.evaluateWithShadow(r, sun, earthRadius, sunRadius);
    };
    require(radiation(0, {-20000, 0, 0}, {}).acceleration.isZero(),
            "SRP must vanish in umbra.");
    relativeNear(radiation(0, {20000, 0, 0}, {}).acceleration,
                 srp.evaluateAtSunPosition({20000, 0, 0}, sun).acceleration, 1e-14);
    checkPartials(radiation, edge, .01);
    rejects([&] { (void)evaluateEclipse({0, 0, 0}, sun, earthRadius, sunRadius); });
    rejects([&] { (void)evaluateEclipse(edge, sun, -1, sunRadius); });

    const auto forces = sumAccelerations({pointMassGravity(398600.4418), radiation});
    fd::filters::CartesianPropagationConfig config;
    config.integrator.absoluteTolerance = config.integrator.relativeTolerance = 1e-12;
    config.integrator.maximumStep = 10;
    Eigen::VectorXd initial(6);
    initial << edge, 0, -4, 0;
    const auto coarse = fd::filters::CartesianPropagator(forces, config)(0, 100, initial);
    config.integrator.maximumStep = 2;
    const fd::filters::CartesianPropagator fine(forces, config);
    const auto reference = fine(0, 100, initial);
    require((coarse.state - reference.state).head<3>().norm() < 1e-6,
            "Eclipse propagation failed maximum-step convergence.");
    require(evaluateEclipse(reference.state.head<3>(), sun, earthRadius, sunRadius).flag ==
                EclipseFlag::Umbra,
            "Propagation did not cross the penumbra into umbra.");
    Eigen::VectorXd plus = initial, minus = initial;
    plus[1] += .01;
    minus[1] -= .01;
    const Eigen::VectorXd numeric = (fine(0, 100, plus).state - fine(0, 100, minus).state) / .02;
    relativeNear(reference.transition.col(1), numeric, 1e-6);
}
} // namespace

int main() {
    try {
        forceLawsAndPartials();
        propagatedPartialsAndNoise();
        eclipseGeometryAndPropagation();
        std::cout
            << "Perturbation laws, analytic partials, composition and RKF45 propagation passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
