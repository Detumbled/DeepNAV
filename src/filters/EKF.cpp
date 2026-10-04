#include "filters/EKF.hpp"
#include "CovarianceValidation.hpp"
#include <Eigen/Cholesky>
#include <cmath>
#include <stdexcept>

namespace fd::filters {

EKF::EKF(StateLayout layout) : layout_(layout) {
    if (layout != StateLayout::Orbit && layout != StateLayout::OrbitClock)
        throw std::invalid_argument("EKF layout must contain 6 or 8 states.");
}

void EKF::setInitialState(const Eigen::VectorXd& state, const Eigen::MatrixXd& covariance,
                          double epoch) {
    if (state.size() != static_cast<int>(layout_))
        throw std::invalid_argument("EKF initial state does not match its layout.");
    detail::validateCovariance(covariance, state.size(), false, "Initial covariance");
    Filter::setInitialState(state, covariance, epoch);
    last_ = {};
}

void EKF::predictTo(double epoch, const PropagationFunction& propagate) {
    if (!initialized_)
        throw std::logic_error("Initialize the EKF before prediction.");
    if (!std::isfinite(epoch) || epoch < epoch_tk_ || !propagate)
        throw std::invalid_argument("EKF prediction requires a model and a finite forward epoch.");
    if (epoch == epoch_tk_)
        return;
    const StatePrediction result = propagate(epoch_tk_, epoch, x_hat_);
    const Eigen::Index n = x_hat_.size();
    if (result.state.size() != n || !result.state.allFinite() || result.transition.rows() != n ||
        result.transition.cols() != n || !result.transition.allFinite())
        throw std::invalid_argument("EKF propagation state/transition is invalid.");
    detail::validateCovariance(result.processCovariance, n, false, "Process covariance");
    const Eigen::MatrixXd covariance = detail::symmetrize(
        result.transition * P_ * result.transition.transpose() + result.processCovariance);
    detail::validateCovariance(covariance, n, false, "Predicted covariance");
    x_hat_ = result.state;
    P_ = covariance;
    epoch_tk_ = epoch;
    last_ = {};
}

InnovationDiagnostics EKF::update(double epoch, const Eigen::VectorXd& observed,
                                  const Eigen::MatrixXd& covariance,
                                  const MeasurementFunction& model, double nisLimit) {
    if (!initialized_)
        throw std::logic_error("Initialize the EKF before an update.");
    if (!std::isfinite(epoch) || epoch != epoch_tk_ || !model || !observed.allFinite())
        throw std::invalid_argument("EKF measurement needs a model and the current epoch.");
    const MeasurementPrediction prediction = model(epoch, x_hat_);
    if (prediction.value.size() != observed.size() || !prediction.value.allFinite())
        throw std::invalid_argument("EKF predicted measurement is invalid.");
    return correct(observed - prediction.value, prediction.jacobian, covariance, nisLimit);
}

void EKF::processBatch(const Eigen::VectorXd& residuals, const Eigen::MatrixXd& jacobian,
                       const Eigen::MatrixXd& covariance) {
    (void)correct(residuals, jacobian, covariance, std::numeric_limits<double>::infinity());
}

InnovationDiagnostics EKF::correct(const Eigen::VectorXd& residuals,
                                   const Eigen::MatrixXd& jacobian,
                                   const Eigen::MatrixXd& covariance, double nisLimit) {
    if (!initialized_)
        throw std::logic_error("Initialize the EKF before an update.");
    const Eigen::Index m = residuals.size(), n = x_hat_.size();
    if (m == 0 || jacobian.rows() != m || jacobian.cols() != n || !residuals.allFinite() ||
        !jacobian.allFinite() || !(nisLimit > 0.0))
        throw std::invalid_argument("EKF measurement dimensions, values or NIS limit are invalid.");
    detail::validateCovariance(covariance, m, true, "Measurement covariance");
    InnovationDiagnostics diagnostics;
    diagnostics.innovation = residuals;
    const Eigen::MatrixXd cross = P_ * jacobian.transpose();
    diagnostics.covariance = detail::symmetrize(jacobian * cross + covariance);
    detail::validateCovariance(diagnostics.covariance, m, true, "Innovation covariance");

    const Eigen::VectorXd scale = diagnostics.covariance.diagonal().array().sqrt().inverse();
    const Eigen::MatrixXd normalized =
        scale.asDiagonal() * diagnostics.covariance * scale.asDiagonal();
    Eigen::LLT<Eigen::MatrixXd> factor(normalized);
    if (factor.info() != Eigen::Success)
        throw std::runtime_error("EKF innovation solve failed.");
    diagnostics.whitenedInnovation = factor.matrixL().solve(scale.asDiagonal() * residuals);
    diagnostics.normalizedInnovationSquared = diagnostics.whitenedInnovation.squaredNorm();
    if (!std::isfinite(diagnostics.normalizedInnovationSquared))
        throw std::runtime_error("EKF innovation statistics overflowed.");
    if (diagnostics.normalizedInnovationSquared > nisLimit) {
        last_ = diagnostics;
        return diagnostics;
    }
    const Eigen::MatrixXd rhs = scale.asDiagonal() * cross.transpose();
    const Eigen::MatrixXd gain = (scale.asDiagonal() * factor.solve(rhs)).transpose();
    const Eigen::VectorXd state = x_hat_ + gain * residuals;
    const Eigen::MatrixXd identityMinusKH = Eigen::MatrixXd::Identity(n, n) - gain * jacobian;
    // Joseph form also remains valid when a future driver deliberately modifies the gain.
    const Eigen::MatrixXd posterior = detail::symmetrize(
        identityMinusKH * P_ * identityMinusKH.transpose() + gain * covariance * gain.transpose());
    if (!state.allFinite())
        throw std::runtime_error("EKF posterior state is nonfinite.");
    detail::validateCovariance(posterior, n, false, "Posterior covariance");
    x_hat_ = state;
    P_ = posterior;
    diagnostics.accepted = true;
    last_ = diagnostics;
    return diagnostics;
}

} // namespace fd::filters
