#include "../benchmarks/cassini_one_way_ekf/ReceptionTrajectory.hpp"
#include "dynamics/LinearStateProvider.hpp"
#include "perturbations/Gravitational.hpp"
#include <iostream>
#include <stdexcept>

int main() {
  using State = Eigen::Matrix<double, 6, 1>;
  const double epoch = 1.5e8;
  const auto force = fd::perturbations::pointMassGravity(398600.);
  const fd::dynamics::LinearStateProvider origin(
      fd::opnav::TdbEpoch{epoch}, {Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()}, "J2000",
      "SSB");
  Eigen::VectorXd initial(8);
  initial << 10000, 2000, 1000, -1, 6, 1, 2e-6, 2e-12;
  const cassini_benchmark::ReceptionTrajectory local(epoch, initial, force,
                                                     origin, 30);
  auto perturbed = initial;
  perturbed.head<3>() += Eigen::Vector3d(.1, -.07, .04);
  perturbed.segment<3>(3) += Eigen::Vector3d(1e-5, -8e-6, 6e-6);
  const auto provider = local.factory()(epoch, perturbed);
  od::RKF45Integrator::Options options;
  options.absoluteTolerance = 1e-13;
  options.relativeTolerance = 1e-13;
  options.maximumStep = 10;
  for (double dt : {-60., -.01, 0., .01}) {
    const auto rate = [&](double t, od::RKF45Integrator::ConstStateRef x,
                          od::RKF45Integrator::StateRef dx) {
      dx.head<3>() = x.tail<3>();
      dx.tail<3>() = force(epoch + t, x.head<3>(), x.tail<3>()).acceleration;
    };
    const auto actual = od::RKF45Integrator(options)
                            .integrate(0, perturbed.head<6>(), dt, rate)
                            .state;
    const auto affine = provider->stateAt(fd::opnav::TdbEpoch{epoch + dt});
    if ((affine.positionKm - actual.head<3>()).norm() > 5e-5 ||
        (affine.velocityKmPerSec - actual.tail<3>()).norm() > 1e-7)
      throw std::runtime_error(
          "Local count trajectory/STM differs from nonlinear propagation");
  }
  std::cout
      << "Forward/backward reception window and affine sensitivity verified.\n";
}
