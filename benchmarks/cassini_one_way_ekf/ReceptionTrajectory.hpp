#pragma once
#include "RKF45Integrator.hpp"
#include "dynamics/EphemerisInterpolator.hpp"
#include "dynamics/StateProvider.hpp"
#include "observations/radiometric/TrajectoryFactory.hpp"
#include "perturbations/ForceModel.hpp"
#include <memory>

namespace cassini_benchmark {
// Short orbit/STM history around an EKF epoch. Only the receiver trajectory is
// needed around the 60 s count; the retarded emitter is the ground station.
class ReceptionTrajectory {
  double epoch_;
  Eigen::Matrix<double, 6, 1> initial_;
  const fd::dynamics::StateProvider &origin_;
  od::EphemerisInterpolator history_;
  const Eigen::Matrix<double, 6, 1> scales_ =
      (Eigen::Matrix<double, 6, 1>() << 1e4, 1e4, 1e4, 10, 10, 10).finished();
  class Affine final : public fd::dynamics::StateProvider {
    const ReceptionTrajectory &history_;
    Eigen::Matrix<double, 6, 1> delta_;

  public:
    Affine(const ReceptionTrajectory &h, const Eigen::VectorXd &x)
        : history_(h), delta_(x.head<6>() - h.initial_) {}
    fd::dynamics::CartesianState stateAt(fd::opnav::TdbEpoch t) const override {
      const auto node =
          history_.history_.interpolate(t.secondsPastJ2000 - history_.epoch_);
      const Eigen::Matrix<double, 6, 6> phi =
          history_.scales_.asDiagonal() *
          Eigen::Map<const Eigen::Matrix<double, 6, 6>>(node.data() + 6) *
          history_.scales_.cwiseInverse().asDiagonal();
      const Eigen::Matrix<double, 6, 1> x =
          history_.scales_.asDiagonal() * node.head<6>() + phi * delta_;
      const auto origin = history_.origin_.stateAt(t);
      return {x.head<3>() + origin.positionKm,
              x.tail<3>() + origin.velocityKmPerSec};
    }
    std::string_view referenceFrame() const noexcept override {
      return "J2000";
    }
    std::string_view referenceOrigin() const noexcept override { return "SSB"; }
  };

public:
  ReceptionTrajectory(double epoch, const Eigen::VectorXd &x,
                      const fd::perturbations::AccelerationFunction &force,
                      const fd::dynamics::StateProvider &origin, double maxStep)
      : epoch_(epoch), initial_(x.head<6>()), origin_(origin) {
    Eigen::VectorXd a = Eigen::VectorXd::Zero(42);
    a.head<6>() = scales_.cwiseInverse().asDiagonal() * initial_;
    Eigen::Map<Eigen::Matrix<double, 6, 6>>(a.data() + 6).setIdentity();
    const auto rate = [&](double dt, od::RKF45Integrator::ConstStateRef s,
                          od::RKF45Integrator::StateRef ds) {
      const Eigen::Matrix<double, 6, 1> physical =
          scales_.asDiagonal() * s.head<6>();
      const auto f = force(epoch + dt, physical.head<3>(), physical.tail<3>());
      Eigen::Matrix<double, 6, 1> dx;
      dx << physical.tail<3>(), f.acceleration;
      ds.head<6>() = scales_.cwiseInverse().asDiagonal() * dx;
      Eigen::Matrix<double, 6, 6> jac = Eigen::Matrix<double, 6, 6>::Zero();
      jac.topRightCorner<3, 3>().setIdentity();
      jac.bottomLeftCorner<3, 3>() = f.positionJacobian;
      jac.bottomRightCorner<3, 3>() = f.velocityJacobian;
      Eigen::Map<Eigen::Matrix<double, 6, 6>>(ds.data() + 6) =
          scales_.cwiseInverse().asDiagonal() * jac * scales_.asDiagonal() *
          Eigen::Map<const Eigen::Matrix<double, 6, 6>>(s.data() + 6);
    };
    od::RKF45Integrator::Options options;
    options.absoluteTolerance = 1e-13;
    options.relativeTolerance = 1e-13;
    options.initialStep = 10;
    options.maximumStep = maxStep;
    for (double end : {-120., 1.}) {
      const auto solution =
          od::RKF45Integrator(options).integrate(0, a, end, rate);
      for (const auto &node : solution.history)
        history_.addNode(node.tdb, node.state, node.derivative);
    }
  }
  fd::observations::radiometric::TrajectoryFactory factory() const {
    return [this](double epoch, const Eigen::VectorXd &x) {
      if (epoch != epoch_)
        throw std::invalid_argument("Reception trajectory epoch mismatch");
      return std::make_unique<Affine>(*this, x);
    };
  }
};
} // namespace cassini_benchmark
