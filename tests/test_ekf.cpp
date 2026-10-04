#include "filters/CartesianPropagator.hpp"
#include "filters/EKF.hpp"
#include "observations/GeometricRadiometricModel.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/SRP.hpp"
#include <Eigen/Cholesky>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace fd::filters;
namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}
void near(const Eigen::MatrixXd& actual, const Eigen::MatrixXd& expected,
          double tolerance = 1e-10) {
    require(actual.rows() == expected.rows() && actual.cols() == expected.cols() &&
                actual.allFinite(),
            "Numerical comparison has invalid dimensions or values.");
    const double error = (actual - expected).norm() / std::max(1.0, expected.norm());
    if (error > tolerance) {
        std::cerr << "Relative error " << error << ", tolerance " << tolerance << '\n';
        throw std::runtime_error("Numerical comparison failed.");
    }
}
template <class Exception = std::invalid_argument, class Function> void rejects(Function action) {
    try {
        action();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error("Expected exception was not thrown.");
}
AccelerationFunction freeMotion() {
    return [](double, const Eigen::Vector3d&, const Eigen::Vector3d&) {
        return AccelerationEvaluation{};
    };
}

void linearReference() {
    CartesianPropagationConfig config;
    config.integrator.maximumStep = 2;
    config.accelerationDiffusion = 1e-8 * Eigen::Matrix3d::Identity();
    const CartesianPropagator propagate(freeMotion(), config);
    Eigen::VectorXd initial(6);
    initial << 10, 20, 30, .1, .2, .3;
    Eigen::MatrixXd prior = Eigen::MatrixXd::Identity(6, 6);
    EKF filter;
    filter.setInitialState(initial, prior, 0);
    filter.predictTo(10, propagate);
    Eigen::MatrixXd phi = Eigen::MatrixXd::Identity(6, 6);
    phi.topRightCorner<3, 3>() = 10 * Eigen::Matrix3d::Identity();
    Eigen::MatrixXd q = Eigen::MatrixXd::Zero(6, 6);
    q.topLeftCorner<3, 3>() = config.accelerationDiffusion * (1000.0 / 3);
    q.topRightCorner<3, 3>() = config.accelerationDiffusion * 50;
    q.bottomLeftCorner<3, 3>() = q.topRightCorner<3, 3>();
    q.bottomRightCorner<3, 3>() = config.accelerationDiffusion * 10;
    const Eigen::VectorXd predicted = phi * initial;
    const Eigen::MatrixXd predictedP = phi * prior * phi.transpose() + q;
    near(filter.state(), predicted);
    near(filter.covariance(), predictedP);

    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(3, 6);
    h.leftCols<3>().setIdentity();
    Eigen::Matrix3d r;
    r << .04, .01, 0, .01, .09, .02, 0, .02, .16;
    Eigen::Vector3d observed = predicted.head<3>() + Eigen::Vector3d(.2, -.1, .3);
    const Eigen::MatrixXd s = h * predictedP * h.transpose() + r;
    const Eigen::MatrixXd gain = s.ldlt().solve(h * predictedP).transpose();
    const Eigen::VectorXd residual = observed - h * predicted;
    const auto diagnostics = filter.update(10, observed, r, [h](double, const Eigen::VectorXd& x) {
        return MeasurementPrediction{h * x, h};
    });
    near(filter.state(), predicted + gain * residual);
    near(filter.covariance(), predictedP - gain * s * gain.transpose());
    require(std::abs(diagnostics.normalizedInnovationSquared -
                     residual.dot(s.ldlt().solve(residual))) < 1e-12,
            "NIS disagrees with analytic reference.");
    require(diagnostics.accepted, "Valid observation was rejected.");

    EKF reordered;
    reordered.setInitialState(predicted, predictedP, 10);
    Eigen::Matrix3d permutation;
    permutation << 0, 0, 1, 1, 0, 0, 0, 1, 0;
    reordered.processBatch(permutation * residual, permutation * h,
                           permutation * r * permutation.transpose());
    near(reordered.state(), filter.state());
    near(reordered.covariance(), filter.covariance());
}

void clockAndCrossCovariance() {
    CartesianPropagationConfig config;
    config.clock = {3e-16 / 86400, 2.25e-26, 1e-30};
    CartesianPropagator propagate(freeMotion(), config);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(8);
    x[6] = 1e-6;
    x[7] = 1e-10;
    Eigen::VectorXd sigmas(8);
    sigmas << 1, 1, 1, .01, .01, .01, 1e-6, 1e-10;
    Eigen::MatrixXd p = sigmas.array().square().matrix().asDiagonal();
    p(0, 6) = p(6, 0) = 2e-7;
    p(3, 7) = p(7, 3) = 2e-13;
    EKF filter(StateLayout::OrbitClock);
    filter.setInitialState(x, p, 0);
    const auto prediction = propagate(0, 120, x);
    const fd::clocks::ClockModel clock(config.clock);
    near(prediction.transition.bottomRightCorner<2, 2>(), clock.transition(120));
    const auto expectedQ = clock.processNoise(120);
    require((prediction.processCovariance.bottomRightCorner<2, 2>() - expectedQ).norm() < 1e-40,
            "Clock diffusion was lost at orbital scales.");
    filter.predictTo(120, propagate);
    near(filter.covariance(), prediction.transition * p * prediction.transition.transpose() +
                                  prediction.processCovariance);
    require(std::abs(filter.state()[6] -
                     (1e-6 + 120e-10 + .5 * config.clock.frequency_drift_per_s * 120 * 120)) <
                1e-20,
            "Clock bias propagation failed.");
    require(filter.covariance()(0, 6) != 0, "Prediction discarded orbit-clock cross covariance.");
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(1, 8);
    h(0, 0) = 1;
    h(0, 6) = 299792.458;
    filter.processBatch(Eigen::VectorXd::Constant(1, .1), h, Eigen::MatrixXd::Constant(1, 1, .01));
    require(filter.state()[6] != prediction.state[6], "Joint update did not estimate clock bias.");

    Eigen::MatrixXd r = Eigen::MatrixXd::Zero(2, 2);
    r(0, 0) = 1e-24;
    r(1, 1) = 1e4;
    Eigen::MatrixXd clockH = Eigen::MatrixXd::Zero(2, 8);
    clockH(0, 6) = 1;
    clockH(1, 0) = 1;
    filter.processBatch(Eigen::Vector2d(1e-12, .1), clockH, r);
    require(std::isfinite(filter.lastInnovation().normalizedInnovationSquared),
            "Mixed-unit innovation solve failed.");
}

void nonlinearPropagation() {
    constexpr double mu = 398600.4418, radius = 20000;
    Eigen::VectorXd x(6);
    x << radius, 0, 0, 0, std::sqrt(mu / radius), 0;
    CartesianPropagationConfig config;
    config.integrator.maximumStep = 2;
    CartesianPropagator propagate(fd::perturbations::pointMassGravity(mu), config);
    const auto prediction = propagate(0, 300, x);
    const double omega = std::sqrt(mu / (radius * radius * radius)), angle = 300 * omega;
    Eigen::VectorXd exact(6);
    exact << radius * std::cos(angle), radius * std::sin(angle), 0,
        -radius * omega * std::sin(angle), radius * omega * std::cos(angle), 0;
    near(prediction.state, exact, 1e-12);
    Eigen::MatrixXd finiteDifference(6, 6);
    for (int j = 0; j < 6; ++j) {
        const double step = j < 3 ? .1 : 1e-4;
        Eigen::VectorXd plus = x, minus = x;
        plus[j] += step;
        minus[j] -= step;
        finiteDifference.col(j) =
            (propagate(0, 300, plus).state - propagate(0, 300, minus).state) / (2 * step);
    }
    near(prediction.transition, finiteDifference, 2e-8);
    const auto split = propagate(120, 300, propagate(0, 120, x).state);
    near(split.state, prediction.state, 1e-12);

    config.accelerationDiffusion = 1e-12 * Eigen::Matrix3d::Identity();
    config.integrator.maximumStep = 4;
    const auto coarse =
        CartesianPropagator(fd::perturbations::pointMassGravity(mu), config)(0, 300, x);
    config.integrator.maximumStep = 1;
    const auto fine =
        CartesianPropagator(fd::perturbations::pointMassGravity(mu), config)(0, 300, x);
    require((fine.processCovariance - coarse.processCovariance).norm() <
                1e-8 * fine.processCovariance.norm(),
            "Integrated process covariance depends excessively on RKF45 maximum step.");
}

void measurementPartials() {
    using namespace fd::observations;
    const fd::dynamics::CartesianState station{{100, -20, 30}, {.1, -.2, .3}};
    for (int n : {6, 8})
        for (auto direction : {LinkDirection::Uplink, LinkDirection::Downlink}) {
            Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
            x.head<6>() << 20000, 4000, 2000, -1, 3, .4;
            if (n == 8) {
                x[6] = 2e-6;
                x[7] = 1e-10;
            }
            const auto model = geometricRadiometricPrediction(x, station, direction, 1e-6, 2e-11);
            for (int j = 0; j < n; ++j) {
                const double step = j < 3 ? .1 : j < 6 ? 1e-5 : j == 6 ? 1e-8 : 1e-10;
                Eigen::VectorXd plus = x, minus = x;
                plus[j] += step;
                minus[j] -= step;
                const Eigen::Vector2d numeric =
                    (geometricRadiometricPrediction(plus, station, direction, 1e-6, 2e-11).value -
                     geometricRadiometricPrediction(minus, station, direction, 1e-6, 2e-11).value) /
                    (2 * step);
                near(numeric, model.jacobian.col(j), 1e-7);
            }
        }
}

void sequentialOrbitRecovery() {
    const CartesianPropagator propagate(fd::perturbations::pointMassGravity(398600.4418));
    const std::array<fd::dynamics::CartesianState, 4> stations{
        {{{6378, 0, 0}, {0, 0, 0}},
         {{0, 6378, 0}, {0, 0, 0}},
         {{0, 0, 6378}, {0, 0, 0}},
         {{-3682, -3682, -3682}, {0, 0, 0}}}};
    const auto model = [&stations](double, const Eigen::VectorXd& x) {
        MeasurementPrediction prediction{Eigen::VectorXd(8), Eigen::MatrixXd(8, x.size())};
        for (int i = 0; i < 4; ++i) {
            const auto pair = fd::observations::geometricRadiometricPrediction(x, stations[i]);
            prediction.value.segment<2>(2 * i) = pair.value;
            prediction.jacobian.middleRows<2>(2 * i) = pair.jacobian;
        }
        return prediction;
    };
    for (int n : {6, 8}) {
        Eigen::VectorXd truth = Eigen::VectorXd::Zero(n);
        truth.head<6>() << 20000, 0, 0, 0, 4, 2;
        Eigen::VectorXd initial = truth;
        initial.head<6>() +=
            (Eigen::Matrix<double, 6, 1>() << .8, -.6, .4, 1e-4, -1e-4, 1e-4).finished();
        Eigen::VectorXd sigma = Eigen::VectorXd::Constant(n, 1e-3);
        sigma.head<3>().setOnes();
        if (n == 8) {
            initial[6] = 5e-7;
            initial[7] = 3e-11;
            sigma[6] = 1e-6;
            sigma[7] = 1e-10;
        }
        EKF filter(n == 6 ? StateLayout::Orbit : StateLayout::OrbitClock);
        filter.setInitialState(initial, sigma.array().square().matrix().asDiagonal(), 0);
        Eigen::MatrixXd r = Eigen::MatrixXd::Zero(8, 8);
        for (int i = 0; i < 4; ++i) {
            r(2 * i, 2 * i) = 1e-4;
            r(2 * i + 1, 2 * i + 1) = 1e-14;
        }
        for (int k = 0; k <= 60; ++k) {
            const double time = k * 60;
            if (k > 0) {
                truth = propagate(time - 60, time, truth).state;
                filter.predictTo(time, propagate);
            }
            // A full hour without measurements exercises forward covariance propagation.
            if (k < 10 || k >= 50)
                (void)filter.update(time, model(time, truth).value, r, model);
            require(filter.state().allFinite() && filter.covariance().allFinite(),
                    "Sequential recovery became nonfinite.");
        }
        require((filter.state() - truth).head<3>().norm() < .01,
                "Noiseless radiometric orbit did not recover after a gap.");
        require((filter.state() - truth).segment<3>(3).norm() < 1e-6,
                "Noiseless velocity did not recover after a gap.");
        if (n == 8)
            require(std::abs(filter.state()[6] - truth[6]) < 1e-8,
                    "Noiseless clock bias did not recover.");
    }
}

void invalidInputsAndGating() {
    const auto model = [](double, const Eigen::VectorXd& x) {
        Eigen::MatrixXd h = Eigen::MatrixXd::Zero(1, x.size());
        h(0, 0) = 1;
        return MeasurementPrediction{h * x, h};
    };
    EKF filter;
    rejects<std::logic_error>([&] { filter.predictTo(1, CartesianPropagator(freeMotion())); });
    rejects([&] { EKF bad(static_cast<StateLayout>(7)); });
    rejects([&] {
        filter.setInitialState(Eigen::VectorXd::Zero(8), Eigen::MatrixXd::Identity(8, 8), 0);
    });
    Eigen::MatrixXd p = Eigen::MatrixXd::Identity(6, 6);
    p(0, 1) = p(1, 0) = 2;
    rejects([&] { filter.setInitialState(Eigen::VectorXd::Zero(6), p, 0); });
    filter.setInitialState(Eigen::VectorXd::Zero(6), Eigen::MatrixXd::Identity(6, 6), 0);
    const auto x0 = filter.state().eval();
    const auto p0 = filter.covariance().eval();
    const Eigen::MatrixXd r = Eigen::MatrixXd::Identity(1, 1);
    const auto gate = filter.update(0, Eigen::VectorXd::Constant(1, 100), r, model, 9);
    require(!gate.accepted && gate.normalizedInnovationSquared == 5000, "NIS gate failed.");
    near(filter.state(), x0);
    near(filter.covariance(), p0);
    rejects([&] { (void)filter.update(1, Eigen::VectorXd::Zero(1), r, model); });
    rejects([&] { (void)filter.update(0, Eigen::VectorXd::Zero(1), -r, model); });
    rejects([&] {
        (void)filter.update(0, Eigen::VectorXd::Zero(1), r, model,
                            std::numeric_limits<double>::quiet_NaN());
    });
    rejects([&] { filter.predictTo(-1, CartesianPropagator(freeMotion())); });
    rejects([&] {
        filter.predictTo(1, [](double, double, const Eigen::VectorXd& x) {
            return StatePrediction{x, Eigen::MatrixXd::Identity(6, 6),
                                   -Eigen::MatrixXd::Identity(6, 6)};
        });
    });
    near(filter.state(), x0);
    near(filter.covariance(), p0);
    require(filter.epoch() == 0, "Rejected input advanced epoch.");
    Eigen::VectorXd sigma(8);
    sigma << 1, 1, 1, 1, 1, 1, 1e-6, 1e-10;
    Eigen::MatrixXd mixed = sigma.array().square().matrix().asDiagonal();
    mixed(6, 7) = mixed(7, 6) = 2e-16;
    EKF eight(StateLayout::OrbitClock);
    rejects([&] { eight.setInitialState(Eigen::VectorXd::Zero(8), mixed, 0); });
    EKF exact;
    exact.setInitialState(Eigen::VectorXd::Zero(6), Eigen::MatrixXd::Zero(6, 6), 0);
    (void)exact.update(0, Eigen::VectorXd::Ones(1), r, model);
    near(exact.state(), Eigen::VectorXd::Zero(6));
}
} // namespace

int main() {
    try {
        linearReference();
        std::cout << "Linear reference passed.\n";
        clockAndCrossCovariance();
        std::cout << "Clock and cross covariance passed.\n";
        nonlinearPropagation();
        std::cout << "Nonlinear propagation passed.\n";
        measurementPartials();
        std::cout << "Measurement partials passed.\n";
        sequentialOrbitRecovery();
        std::cout << "Sequential orbit recovery passed.\n";
        invalidInputsAndGating();
        std::cout << "EKF reference, clock, Jacobian, covariance and validation checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
