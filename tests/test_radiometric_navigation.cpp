#include "dynamics/IntegratedStateProvider.hpp"
#include "dynamics/LinearStateProvider.hpp"
#include "observations/radiometric/TwoWayMeasurement.hpp"
#include "perturbations/Gravitational.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace radio = fd::observations::radiometric;
namespace {
class CountingProvider final : public fd::dynamics::StateProvider {
  public:
    explicit CountingProvider(double range) : range_(range) {}
    fd::dynamics::CartesianState stateAt(fd::opnav::TdbEpoch) const override {
        ++calls;
        return {{range_, 0, 0}, {0, 0, 0}};
    }
    std::string_view referenceFrame() const noexcept override { return "J2000"; }
    std::string_view referenceOrigin() const noexcept override { return "SSB"; }
    mutable int calls{};

  private:
    double range_;
};
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace
int main() {
    try {
        CountingProvider tx(3e6), rx(0);
        const auto seed = fd::opnav::LightTimeSolver{}.solve(
            fd::opnav::TdbEpoch{100}, rx.stateAt(fd::opnav::TdbEpoch{100}), tx);
        const int geometricCalls = tx.calls;
        tx.calls = rx.calls = 0;
        const auto geometric = radio::libraryReceptionSolver()(100, tx, rx);
        require(rx.calls == 1 && tx.calls == geometricCalls,
                "Geometric endpoint states were queried again");
        require(geometric.emissionTdb == seed.emissionEpoch.secondsPastJ2000,
                "Seed emission changed");
        tx.calls = rx.calls = 0;
        const auto delayed = radio::libraryReceptionSolver(
            {}, [](const auto &, const auto &, const auto &) { return 1e-5; })(100, tx, rx);
        require(rx.calls == 1 && delayed.residualSeconds <= 1e-9,
                "Corrected receiver was not fixed");
        const double mu = 398600.4418, radius = 20000, speed = std::sqrt(mu / radius),
                     epoch = 1.5e8;
        Eigen::Matrix<double, 6, 1> state;
        state << radius, 0, 0, 0, speed, 0;
        od::RKF45Integrator::Options options;
        options.maximumStep = 30;
        options.relativeTolerance = 1e-12;
        options.absoluteTolerance = 1e-11;
        const fd::dynamics::IntegratedStateProvider orbit(
            epoch, state, epoch + 600, fd::perturbations::pointMassGravity(mu), options);
        const double t = 317.3, angle = std::sqrt(mu / std::pow(radius, 3)) * t;
        const Eigen::Vector3d expected(radius * std::cos(angle), radius * std::sin(angle), 0);
        require((orbit.stateAt(fd::opnav::TdbEpoch{epoch + t}).positionKm - expected).norm() < 1e-5,
                "Dense trajectory failed circular-orbit accuracy check");
        const fd::dynamics::LinearStateProvider station(fd::opnav::TdbEpoch{0},
                                                        {{0, 0, 0}, {0, 0, 0}}, "J2000", "SSB");
        const radio::TrajectoryFactory factory = [](double t, const Eigen::VectorXd &x) {
            return std::make_unique<fd::dynamics::LinearStateProvider>(
                fd::opnav::TdbEpoch{t}, fd::dynamics::CartesianState{x.head<3>(), x.tail<3>()},
                "J2000", "SSB");
        };
        Eigen::VectorXd steps(6);
        steps << .1, .1, .1, 1e-5, 1e-5, 1e-5;
        const double ratio = 1.25, f = 2e9;
        const radio::TwoWayMeasurement model(
            radio::TwoWayLink(radio::libraryReceptionSolver(), {ratio, 0}), station, 100, 60,
            {0, f, 0}, ratio * f, {}, factory, steps);
        Eigen::VectorXd x = Eigen::VectorXd::Zero(6);
        x[0] = 3e6;
        const auto prediction = model(0, x);
        require(std::abs(prediction.jacobian(0, 0) - 1) < 1e-7, "Two-way range partial failed");
        require(std::abs(prediction.jacobian(1, 3) +
                         2 * ratio * f / radio::speedOfLightKmPerSecond) < .02,
                "Coherent Doppler velocity partial failed");
        std::cout << "State reuse, corrected convergence, dense trajectory and two-way Jacobian "
                     "passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
