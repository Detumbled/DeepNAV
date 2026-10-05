#pragma once

#include <Eigen/Cholesky>
#include <cmath>
#include <random>
#include <stdexcept>

namespace fd::montecarlo {

// Draw the complete correlated prior, with normalization for orbit/clock units.
inline Eigen::VectorXd drawError(const Eigen::MatrixXd& covariance, std::mt19937_64& generator) {
    const Eigen::VectorXd sigma = covariance.diagonal().array().sqrt();
    if (covariance.rows() == 0 || covariance.rows() != covariance.cols() ||
        !covariance.allFinite() || !sigma.allFinite() || (sigma.array() <= 0).any())
        throw std::invalid_argument("Sampling requires a positive-definite covariance.");
    const Eigen::MatrixXd correlation =
        sigma.cwiseInverse().asDiagonal() * covariance * sigma.cwiseInverse().asDiagonal();
    if (!correlation.allFinite() || !correlation.isApprox(correlation.transpose(), 1e-12))
        throw std::invalid_argument("Covariance must be finite and symmetric after scaling.");
    Eigen::LLT<Eigen::MatrixXd> factor(correlation);
    if (factor.info() != Eigen::Success)
        throw std::runtime_error("Prior covariance factorization failed.");
    std::normal_distribution<double> normal;
    Eigen::VectorXd z(sigma.size());
    for (auto& value : z)
        value = normal(generator);
    return sigma.asDiagonal() * (factor.matrixL() * z);
}

// Full covariance quadratic form, not a sum of marginal error/sigma ratios.
inline double nees(const Eigen::VectorXd& error, const Eigen::MatrixXd& covariance) {
    if (error.size() == 0 || covariance.rows() != error.size() ||
        covariance.cols() != error.size() || !error.allFinite() || !covariance.allFinite() ||
        (covariance.diagonal().array() <= 0).any())
        throw std::invalid_argument(
            "NEES requires matching error and positive-definite covariance.");
    const Eigen::VectorXd scale = covariance.diagonal().array().sqrt().inverse();
    const Eigen::MatrixXd correlation = scale.asDiagonal() * covariance * scale.asDiagonal();
    if (!correlation.allFinite() || !correlation.isApprox(correlation.transpose(), 1e-12))
        throw std::invalid_argument("Covariance must be finite and symmetric after scaling.");
    Eigen::LLT<Eigen::MatrixXd> factor(correlation);
    if (factor.info() != Eigen::Success)
        throw std::runtime_error("NEES covariance factorization failed.");
    const Eigen::VectorXd white = factor.matrixL().solve(scale.asDiagonal() * error);
    const double value = white.squaredNorm();
    if (!std::isfinite(value))
        throw std::runtime_error("NEES overflowed.");
    return value;
}

} // namespace fd::montecarlo
