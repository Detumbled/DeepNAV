#pragma once
#include "filters/BatchLeastSquaresDriver.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace fd::simulation {
struct BatchPrediction {
  Eigen::VectorXd measurements;
  std::vector<Eigen::VectorXd> states;
};
struct BatchBaseline {
  Eigen::VectorXd state;
  Eigen::MatrixXd covariance;
  BatchPrediction prediction;
  std::vector<Eigen::MatrixXd> transitions;
  bool converged{};
  std::size_t iterations{};
  double undampedCorrectionNorm{};
};
using BatchEvaluator = std::function<BatchPrediction(const Eigen::VectorXd &)>;
// Reuse WLS and its iterative driver, in prior-standardized coordinates. All
// observations refer to one initial nine-state vector and a deterministic arc.
inline BatchBaseline solveOneWayBatch(const Eigen::VectorXd &initial,
                                      const Eigen::VectorXd &scales,
                                      const Eigen::VectorXd &observed,
                                      const Eigen::VectorXd &variances,
                                      const BatchEvaluator &evaluate,
                                      const std::filesystem::path &output,
                                      int maxIterations, double arcSeconds) {
  filters::WLS wls;
  wls.setInitialState(Eigen::VectorXd::Zero(9), Eigen::MatrixXd::Identity(9, 9),
                      0);
  Eigen::VectorXd stencil(9);
  stencil << 500, 500, 500, 10, 10, 10, 1e7, 10, 1000;
  // A short-arc SRP perturbation must exceed barycentric rounding. The force
  // is linear in CR; retain a wide stencil even on long arcs to resolve
  // its small COUNTED Doppler sensitivity.
  stencil[8] = std::max(
      1000., 200 * std::pow(6000 / std::max(6000., arcSeconds + 5000), 2));
  const auto physical = [&](const Eigen::VectorXd &z) -> Eigen::VectorXd {
    return initial + scales.cwiseProduct(z);
  };
  const auto linearizeRaw = [&](const Eigen::VectorXd &z) {
    const auto center = evaluate(physical(z));
    filters::BatchLeastSquaresDriver::LinearizedProblem p;
    p.residuals = observed - center.measurements;
    p.measurementVariances = variances;
    p.designMatrix.resize(observed.size(), 9);
    for (int j = 0; j < 9; ++j) {
      auto plus = z.eval(), minus = z.eval();
      plus[j] += stencil[j];
      minus[j] -= stencil[j];
      p.designMatrix.col(j) = (evaluate(physical(plus)).measurements -
                               evaluate(physical(minus)).measurements) /
                              (2 * stencil[j]);
    }
    return p;
  };
  filters::BatchLeastSquaresDriver::ConvergenceCriteria criteria;
  criteria.maxIterations = maxIterations;
  criteria.minimumIterations = 2;
  criteria.correctionComponentTolerances = Eigen::VectorXd::Constant(9, 1e-4);
  const auto objective = [&](const Eigen::VectorXd &z,
                             const Eigen::VectorXd &residual) {
    return (residual.array().square() / variances.array()).sum() +
           z.squaredNorm();
  };
  double damping = 0;
  const auto linearize = [&](const Eigen::VectorXd &z) {
    auto raw = linearizeRaw(z);
    const double current = objective(z, raw.residuals);
    for (int attempt = 0; attempt < 14; ++attempt) {
      auto candidate = raw;
      if (damping > 0) {
        candidate.residuals.conservativeResize(observed.size() + 9);
        candidate.residuals.tail<9>().setZero();
        candidate.designMatrix.conservativeResize(observed.size() + 9, 9);
        candidate.designMatrix.bottomRows<9>() =
            std::sqrt(damping) * Eigen::MatrixXd::Identity(9, 9);
        candidate.measurementVariances.conservativeResize(observed.size() + 9);
        candidate.measurementVariances.tail<9>().setOnes();
      }
      auto probe = wls;
      probe.processBatchDiagonal(candidate.residuals, candidate.designMatrix,
                                 candidate.measurementVariances);
      const auto predicted = evaluate(physical(probe.state()));
      if (objective(probe.state(), observed - predicted.measurements) <=
          current + 1e-8)
        return candidate;
      damping = damping == 0 ? 1 : damping * 10;
    }
    throw std::runtime_error(
        "WLS could not find a decreasing regularized objective");
  };
  filters::BatchLeastSquaresDriver driver(wls, linearize, criteria);
  driver.setConvergenceFunction([&](const auto &it, const auto &c) {
    const auto before = evaluate(physical(it.stateBefore)),
               after = evaluate(physical(it.stateAfter));
    const double a = objective(it.stateBefore, observed - before.measurements);
    const double b = objective(it.stateAfter, observed - after.measurements);
    const bool smallCorrection = it.correction.cwiseAbs().maxCoeff() < 1e-4;
    const bool stationary =
        std::abs(a - b) / std::max(1., a) < 1e-7 && it.correction.norm() < .01;
    damping *= .1;
    return it.iteration >= c.minimumIterations &&
           (smallCorrection || stationary);
  });
  std::ofstream log(output / "wls_iterations.csv");
  log.exceptions(std::ios::badbit | std::ios::failbit);
  log << std::setprecision(17)
      << "iteration,weighted_prefit_rms,nonlinear_postfit_rms,normalized_"
         "correction_norm,"
         "regularized_objective,damping";
  for (int j = 0; j < 9; ++j)
    log << ",correction_" << j;
  log << '\n';
  const auto result = driver.solve([&](const auto &it) {
    const auto predicted = evaluate(physical(it.stateAfter));
    const Eigen::VectorXd residual = observed - predicted.measurements;
    const double postRms =
        std::sqrt((residual.array().square() / variances.array()).mean());
    log << it.iteration << ',' << it.weightedRms << ',' << postRms << ','
        << it.correctionNorm << ',' << objective(it.stateAfter, residual) << ','
        << damping;
    for (int j = 0; j < 9; ++j)
      log << ',' << it.correction[j];
    log << '\n';
    std::cout << "WLS iteration " << it.iteration
              << ": nonlinear weighted RMS=" << postRms
              << ", normalized correction=" << it.correctionNorm
              << ", damping=" << damping << std::endl;
  });
  BatchBaseline out;
  out.state = physical(wls.state());
  // A small damped step alone cannot prove convergence. The final raw
  // (undamped) correction below independently certifies stationarity.
  out.converged = result.converged();
  out.iterations = result.iterationCount();
  out.prediction = evaluate(out.state);
  out.transitions.assign(out.prediction.states.size(),
                         Eigen::MatrixXd::Zero(9, 9));
  // Recompute final undamped information and trajectory derivatives at the
  // FINAL nonlinear solution, rather than using a previous iteration's H.
  const auto finalProblem = linearizeRaw(wls.state());
  auto covarianceSolver = wls;
  covarianceSolver.processBatchDiagonal(finalProblem.residuals,
                                        finalProblem.designMatrix, variances);
  out.covariance =
      scales.asDiagonal() * covarianceSolver.covariance() * scales.asDiagonal();
  out.undampedCorrectionNorm = covarianceSolver.lastStateCorrection().norm();
  out.converged = out.converged && out.undampedCorrectionNorm < .01;
  for (int j = 0; j < 9; ++j) {
    auto plus = out.state.eval(), minus = out.state.eval();
    const double step = stencil[j] * scales[j];
    plus[j] += step;
    minus[j] -= step;
    const auto a = evaluate(plus), b = evaluate(minus);
    for (std::size_t k = 0; k < out.transitions.size(); ++k)
      out.transitions[k].col(j) = (a.states[k] - b.states[k]) / (2 * step);
  }
  return out;
}
} // namespace fd::simulation
