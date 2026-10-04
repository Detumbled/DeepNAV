#include "filters/CartesianPropagator.hpp"
#include "CovarianceValidation.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::filters {
namespace {
using State6 = Eigen::Matrix<double, 6, 1>;
using Matrix6 = Eigen::Matrix<double, 6, 6>;
} // namespace

CartesianPropagator::CartesianPropagator(AccelerationFunction acceleration,
                                         CartesianPropagationConfig config)
    : acceleration_(std::move(acceleration)), config_(std::move(config)), clock_(config_.clock),
      integrator_(config_.integrator) {
    if (!acceleration_ || !config_.stateScales.allFinite() ||
        (config_.stateScales.array() <= 0).any())
        throw std::invalid_argument(
            "Cartesian propagation requires a force model and positive finite state scales.");
    detail::validateCovariance(config_.accelerationDiffusion, 3, false, "Acceleration diffusion");
}

StatePrediction CartesianPropagator::operator()(double fromEpoch, double toEpoch,
                                                const Eigen::VectorXd& initial) const {
    const double duration = toEpoch - fromEpoch;
    const Eigen::Index n = initial.size();
    if (!std::isfinite(fromEpoch) || !std::isfinite(toEpoch) || !std::isfinite(duration) ||
        duration < 0.0 || (n != 6 && n != 8) || !initial.allFinite())
        throw std::invalid_argument(
            "Cartesian propagation requires finite 6/8-state input and forward time.");
    StatePrediction result{initial, Eigen::MatrixXd::Identity(n, n), Eigen::MatrixXd::Zero(n, n)};
    if (duration == 0.0)
        return result;
    const bool withNoise = (config_.accelerationDiffusion.array() != 0.0).any();
    const State6 inverseScale = config_.stateScales.cwiseInverse();
    State6 noiseScale = State6::Ones();
    if (withNoise) {
        const double velocityScale =
            std::sqrt(config_.accelerationDiffusion.diagonal().maxCoeff()) * std::sqrt(duration);
        noiseScale.head<3>().setConstant(velocityScale * duration);
        noiseScale.tail<3>().setConstant(velocityScale);
        if (!noiseScale.allFinite() || (noiseScale.array() <= 0).any())
            throw std::invalid_argument("Process covariance scales overflowed or underflowed.");
    }
    const State6 inverseNoiseScale = noiseScale.cwiseInverse();
    Eigen::VectorXd augmented = Eigen::VectorXd::Zero(withNoise ? 78 : 42);
    augmented.head<6>() = inverseScale.asDiagonal() * initial.head<6>();
    Eigen::Map<Matrix6>(augmented.data() + 6).setIdentity();
    const auto dynamics = [&](double elapsed, od::RKF45Integrator::ConstStateRef state,
                              od::RKF45Integrator::StateRef rate) {
        const State6 physical = config_.stateScales.asDiagonal() * state.head<6>();
        const auto force =
            acceleration_(fromEpoch + elapsed, physical.head<3>(), physical.tail<3>());
        if (!force.acceleration.allFinite() || !force.positionJacobian.allFinite() ||
            !force.velocityJacobian.allFinite())
            throw std::invalid_argument("Cartesian force/partials must be finite.");
        State6 physicalRate;
        physicalRate << physical.tail<3>(), force.acceleration;
        rate.head<6>() = inverseScale.asDiagonal() * physicalRate;
        Matrix6 jacobian = Matrix6::Zero();
        jacobian.topRightCorner<3, 3>().setIdentity();
        jacobian.bottomLeftCorner<3, 3>() = force.positionJacobian;
        jacobian.bottomRightCorner<3, 3>() = force.velocityJacobian;
        Eigen::Map<Matrix6>(rate.data() + 6) = inverseScale.asDiagonal() * jacobian *
                                               config_.stateScales.asDiagonal() *
                                               Eigen::Map<const Matrix6>(state.data() + 6);
        if (withNoise) {
            const Matrix6 a = inverseNoiseScale.asDiagonal() * jacobian * noiseScale.asDiagonal();
            const Eigen::Map<const Matrix6> q(state.data() + 42);
            Eigen::Map<Matrix6> qRate(rate.data() + 42);
            // Integrate the discrete noise integral with the same physical dynamics as the STM.
            qRate = a * q + q * a.transpose();
            qRate.bottomRightCorner<3, 3>() += inverseNoiseScale.tail<3>().asDiagonal() *
                                               config_.accelerationDiffusion *
                                               inverseNoiseScale.tail<3>().asDiagonal();
        }
    };
    // Integrate elapsed time to avoid adding tiny adaptive steps to large TDB epochs.
    const auto integrated = integrator_.integrate(0.0, augmented, duration, dynamics);
    result.state.head<6>() = config_.stateScales.asDiagonal() * integrated.state.head<6>();
    result.transition.topLeftCorner<6, 6>() =
        config_.stateScales.asDiagonal() * Eigen::Map<const Matrix6>(integrated.state.data() + 6) *
        inverseScale.asDiagonal();
    if (withNoise)
        result.processCovariance.topLeftCorner<6, 6>() = detail::symmetrize(
            noiseScale.asDiagonal() * Eigen::Map<const Matrix6>(integrated.state.data() + 42) *
            noiseScale.asDiagonal());
    if (n == 8) {
        const auto clock = clock_.propagate({initial[6], initial[7]}, duration);
        result.state[6] = clock.bias_s;
        result.state[7] = clock.fractional_frequency;
        result.transition.bottomRightCorner<2, 2>() = clock_.transition(duration);
        result.processCovariance.bottomRightCorner<2, 2>() = clock_.processNoise(duration);
    }
    detail::validateCovariance(result.processCovariance, n, false, "Integrated process covariance");
    return result;
}

} // namespace fd::filters
