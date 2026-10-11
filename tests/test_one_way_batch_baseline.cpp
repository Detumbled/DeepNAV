#include "filters/OneWayBatchBaseline.hpp"
#include <Eigen/Dense>
#include <stdexcept>

int main() {
  Eigen::VectorXd initial(9), scales(9), truth(9);
  initial << 10000, -8000, 5000, 1, -.8, .5, 2e-6, 2e-11, 1.2;
  scales << 2, 2, 2, 1e-4, 1e-4, 1e-4, 1e-5, 1e-10, .5;
  truth << .5, -.4, .3, .2, -.1, .6, .25, -.2, .7;
  Eigen::MatrixXd design = Eigen::MatrixXd::Zero(18, 9);
  design.topRows<9>().setIdentity();
  design.bottomRows<9>() = 2 * Eigen::MatrixXd::Identity(9, 9);
  design(0, 6) = 3; // coupled orbit/clock parameter
  design(8, 2) = -.5;
  const Eigen::VectorXd variances = Eigen::VectorXd::Constant(18, .25);
  const Eigen::VectorXd observed = design * truth;
  const auto evaluate = [&](const Eigen::VectorXd &x) {
    fd::simulation::BatchPrediction p;
    p.measurements = design * (x - initial).cwiseQuotient(scales);
    p.states = {x, 2 * x};
    return p;
  };
  const auto output =
      std::filesystem::current_path() / "Output_batch_baseline_test";
  std::filesystem::create_directories(output);
  const auto actual = fd::simulation::solveOneWayBatch(
      initial, scales, observed, variances, evaluate, output, 6, 1209600);
  const Eigen::MatrixXd information =
      Eigen::MatrixXd::Identity(9, 9) + 4 * design.transpose() * design;
  const Eigen::VectorXd expected =
      initial + scales.cwiseProduct(information.ldlt().solve(
                    4 * design.transpose() * observed));
  const auto require = [](bool okay, const char *message) {
    if (!okay)
      throw std::runtime_error(message);
  };
  require(actual.converged, "Linear batch failed to converge");
  require(((actual.state - expected).cwiseQuotient(scales)).norm() < 1e-7,
          "Batch state differs from analytic MAP solution");
  const Eigen::MatrixXd normalizedCovariance =
      scales.cwiseInverse().asDiagonal() * actual.covariance *
      scales.cwiseInverse().asDiagonal();
  require((normalizedCovariance - information.inverse()).norm() < 1e-7,
          "Batch covariance differs from analytic posterior");
  require(
      (actual.transitions[0] - Eigen::MatrixXd::Identity(9, 9)).norm() < 1e-6 &&
          (actual.transitions[1] - 2 * Eigen::MatrixXd::Identity(9, 9)).norm() <
              1e-6,
      "Full-arc state/covariance sensitivities are inconsistent");
  // A numerical starting point must not silently become the prior mean.
  for (double start : {-100., 10.}) {
    fd::filters::WLS filter;
    filter.setInitialState(Eigen::VectorXd::Constant(1,start),Eigen::MatrixXd::Constant(1,1,4.),0);
    filter.setAprioriMean(Eigen::VectorXd::Zero(1));
    filter.processBatchDiagonal(Eigen::VectorXd::Constant(1,2.-start),Eigen::MatrixXd::Identity(1,1),Eigen::VectorXd::Ones(1));
    require(std::abs(filter.state()[0]-1.6)<1e-12,"Changing the starting guess changed the MAP estimate");
    require(std::abs(filter.covariance()(0,0)-.8)<1e-12,"Independent prior mean changed the covariance");
    filter.processBatchDiagonal(Eigen::VectorXd::Constant(1,2.-filter.state()[0]),Eigen::MatrixXd::Identity(1,1),Eigen::VectorXd::Ones(1));
    require(std::abs(filter.lastStateCorrection()[0])<1e-12,"Prior mean was lost on relinearization");
  }
  std::filesystem::remove_all(output);
}
