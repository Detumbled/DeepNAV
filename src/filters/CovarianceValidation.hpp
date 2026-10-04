#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <stdexcept>
#include <string>

namespace fd::filters::detail {

inline void validateCovariance(const Eigen::MatrixXd& matrix, Eigen::Index size,
                               bool positiveDefinite, const char* name) {
    const auto fail = [name] {
        throw std::invalid_argument(
            std::string(name) +
            " must be finite, symmetric and positive (semi)definite with matching dimensions.");
    };
    if (size <= 0 || matrix.rows() != size || matrix.cols() != size || !matrix.allFinite())
        fail();
    Eigen::VectorXd scale(size);
    for (Eigen::Index i = 0; i < size; ++i) {
        if (matrix(i, i) < 0.0 || (positiveDefinite && matrix(i, i) == 0.0))
            fail();
        if (matrix(i, i) == 0.0 &&
            ((matrix.row(i).array() != 0.0).any() || (matrix.col(i).array() != 0.0).any()))
            fail();
        scale[i] = matrix(i, i) > 0.0 ? 1.0 / std::sqrt(matrix(i, i)) : 1.0;
    }
    // Normalize units before testing definiteness: clock and orbit scales differ greatly.
    const Eigen::MatrixXd normalized = scale.asDiagonal() * matrix * scale.asDiagonal();
    if (!normalized.allFinite() || !normalized.isApprox(normalized.transpose(), 1e-12))
        fail();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(normalized);
    if (eigen.info() != Eigen::Success || eigen.eigenvalues().minCoeff() < -1e-12)
        fail();
    if (positiveDefinite) {
        Eigen::LLT<Eigen::MatrixXd> factor(normalized);
        if (factor.info() != Eigen::Success)
            fail();
    }
}

inline Eigen::MatrixXd symmetrize(const Eigen::MatrixXd& matrix) {
    return 0.5 * (matrix + matrix.transpose());
}

} // namespace fd::filters::detail
