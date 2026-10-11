#include "dynamics/IntegratedStateProvider.hpp"
#include "dynamics/LinearStateProvider.hpp"
#include "filters/CartesianPropagator.hpp"
#include "filters/EKF.hpp"
#include "observations/radiometric/OneWayMeasurement.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace radio = fd::observations::radiometric;
namespace {
void near(double a, double b, double tol, const char *label) {
    if (!std::isfinite(a) || std::abs(a - b) > tol)
        throw std::runtime_error(std::string(label) + ": " + std::to_string(a) + " vs " +
                                 std::to_string(b));
}
auto linear(double r, double v = 0, double epoch = 0) {
    return fd::dynamics::LinearStateProvider(fd::opnav::TdbEpoch{epoch}, {{r, 0, 0}, {v, 0, 0}},
                                             "J2000", "SSB");
}
void physics() {
    const double f = 7.175e9, c = radio::speedOfLightKmPerSecond;
    const radio::OneWayLink link(radio::libraryReceptionSolver());
    for (double epoch : {0., 1.5e8}) {
        const auto tx = linear(0, 0, epoch), rx = linear(3e6, 10, epoch);
        const radio::EndpointClock ideal{epoch, {}, {}},
            onboard{epoch, [=](double t) { return 2e-6 + 2e-10 * (t - epoch); }, {}};
        auto a = link.count(100, 60, tx, rx, ideal, ideal, {epoch, f, 0}, f);
        near(a.residualHz, -f * 10 / c, .002, "Receding uplink phase Doppler");
        a = link.count(onboard.reading(epoch + 100), 60, tx, rx, ideal, onboard, {epoch, f, 0}, f);
        near(a.residualHz, f * ((1 - 10 / c) / (1 + 2e-10) - 1), .002,
             "Clock frequency sign and count interval");
        near(a.receiveCoordinateIntervalSeconds, 60 / (1 + 2e-10), 1e-11,
             "Count in onboard seconds");
        const auto stationary = linear(3e6, 0, epoch);
        const radio::EndpointClock bias{epoch, [](double) { return 1e-6; }, {}};
        a = link.count(bias.reading(epoch + 100), 60, tx, stationary, ideal, bias, {epoch, f, 0},
                       f);
        near(a.pseudorangeKm, 3e6 + c * 1e-6, 1e-8,
             "One microsecond one-way pseudorange = 299.792458 m");
        near(a.residualHz, 0, 1e-8, "Constant receiver bias does not shift stationary carrier");
        const auto b = link.count(100, 60, tx, stationary, ideal, ideal, {epoch, f, 2}, f);
        near(b.residualHz, 2 * (70 - 3e6 / c), 1e-5, "Ramp at retarded ground transmission");
        const radio::EndpointClock ground{
            epoch, [=](double t) { return 1e-6 + 1e-10 * (t - epoch); }, {}};
        const auto d = link.count(100, 60, tx, stationary, ground, ideal, {epoch, f, 0}, f);
        near(d.residualHz, f * 1e-10, 1e-6, "Ground frequency sign");
        const radio::EndpointClock offset{epoch, {}, [](double) { return 123.; }};
        near(offset.coordinateEpoch(offset.reading(epoch + 100)), epoch + 100, 1e-8,
             "Nominal reference inversion");
        const auto movingGround = linear(0, 10, epoch);
        const auto stationMotion =
            link.count(100, 60, movingGround, stationary, ideal, ideal, {epoch, f, 0}, f);
        near(stationMotion.residualHz, f * (10 / (c - 10)), .002,
             "Retarded moving transmitter Doppler");
    }
    const auto tx = linear(0), rx = linear(3e6);
    const radio::EndpointClock ideal{};
    const radio::OneWayLink dispersive(
        radio::libraryReceptionSolver(
            {}, [](const auto &, const auto &, const auto &) { return 1e-6; }),
        radio::libraryReceptionSolver());
    const auto a = dispersive.count(100, 60, tx, rx, ideal, ideal, {0, f, 0}, f);
    near(a.pseudorangeKm, 3e6 + c * 1e-6, 1e-8, "Separate group correction");
    near(a.end.correctionSeconds, 0, 0, "Separate carrier phase correction");
    const double drift = 1e-7;
    const radio::EndpointClock drifting{0, [=](double t) { return .5 * drift * t * t; }, {}};
    const auto driftCount =
        link.count(drifting.reading(100), 60, tx, rx, ideal, drifting, {0, f, 0}, f);
    const double rate = 1 + drift * 100;
    const double dt = 120 / (rate + std::sqrt(rate * rate - 120 * drift));
    near(driftCount.receiveCoordinateIntervalSeconds, dt, 1e-11,
         "Quadratic clock drift interval inversion");
    near(driftCount.residualHz, f * (dt / 60 - 1), 1e-5, "Quadratic clock drift phase count");
    const radio::EndpointClock nominalRate{0, {}, [](double t) { return .001 + 1e-8 * t; }};
    const auto nominalCount =
        link.count(nominalRate.reading(100), 60, tx, rx, ideal, nominalRate, {0, f, 0}, f);
    near(nominalCount.residualHz, f * 1e-8 / (1 - 1e-8), 1e-6,
         "Time-dependent nominal reference offset");
    bool rejected = false;
    try {
        (void)link.count(100, 0, tx, rx, ideal, ideal, {0, f, 0}, f);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("Zero count accepted");
}
// Analytic accelerated motion isolates weak clock-bias Doppler sensitivity at
// a deep-space distance and a large absolute TDB epoch.
class AcceleratedLine final : public fd::dynamics::StateProvider {
public:
    AcceleratedLine(double epoch, const Eigen::VectorXd &x) : epoch_(epoch), x_(x) {}
    fd::dynamics::CartesianState stateAt(fd::opnav::TdbEpoch t) const override {
        const double dt = t.secondsPastJ2000 - epoch_;
        fd::dynamics::CartesianState s{x_.head<3>(), x_.segment<3>(3)};
        s.positionKm += s.velocityKmPerSec * dt;
        s.positionKm[0] += .5e-6 * dt * dt;
        s.velocityKmPerSec[0] += 1e-6 * dt;
        return s;
    }
    std::string_view referenceFrame() const noexcept override { return "J2000"; }
    std::string_view referenceOrigin() const noexcept override { return "SSB"; }
private:
    double epoch_;
    Eigen::VectorXd x_;
};
void weakClockPartialAtDeepSpaceEpoch() {
    const double epoch = 1.5e8, f = 7.175e9, c = radio::speedOfLightKmPerSecond;
    const auto tx = linear(0, 0, epoch);
    const radio::TrajectoryFactory factory = [](double t, const Eigen::VectorXd &x) {
        return std::make_unique<AcceleratedLine>(t, x);
    };
    Eigen::VectorXd x = Eigen::VectorXd::Zero(9), steps(9);
    x[0] = 1e9; x[3] = 5; x[8] = 1.2;
    steps << 10, 10, 10, 1e-3, 1e-3, 1e-3, 1, 1e-9, 10;
    for (double clockStep : {5., 10., 20.}) {
        steps[6] = clockStep;
        const radio::OneWayMeasurement model(
            radio::OneWayLink(radio::libraryReceptionSolver()), tx, 260, 60, epoch,
            {epoch, f, 0}, f, {epoch, {}, {}}, {}, 0, factory, steps);
        const auto p = model(epoch, x);
        near(p.jacobian(0, 6), c - 5 - 260e-6, .001,
             "Accelerated deep-space range/clock derivative");
        near(p.jacobian(1, 6), f * 1e-6 / c, .0003,
             "Weak Doppler/clock derivative at large epoch");
        near(p.jacobian(1, 0), 0, 1e-5,
             "Constant-acceleration Doppler independent of initial range");
    }
}
void measurementAndNineStates() {
    const auto tx = linear(0);
    const double f = 2e9, c = radio::speedOfLightKmPerSecond;
    const radio::TrajectoryFactory factory = [](double t, const Eigen::VectorXd &x) {
        return std::make_unique<fd::dynamics::LinearStateProvider>(
            fd::opnav::TdbEpoch{t}, fd::dynamics::CartesianState{x.head<3>(), x.segment<3>(3)},
            "J2000", "SSB");
    };
    Eigen::VectorXd steps(9);
    steps << .1, .1, .1, 1e-5, 1e-5, 1e-5, 1e-6, 1e-9, .1;
    const radio::OneWayMeasurement model(radio::OneWayLink(radio::libraryReceptionSolver()), tx,
                                         100, 60, 0, {0, f, 0}, f, {}, {}, 0, factory, steps);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(9);
    x[0] = 3e6;
    x[8] = 1.2;
    const auto p = model(0, x);
    near(p.jacobian(0, 0), 1, 1e-7, "Range orbit partial");
    near(p.jacobian(0, 6), c, .002, "Range clock bias partial");
    near(p.jacobian(0, 7), 100 * c, 1, "Range clock frequency partial");
    near(p.jacobian(1, 7), -f, 10, "Carrier receiver frequency partial");
    near(p.jacobian(1, 3), -f / c, .1, "One-way velocity partial");
    near(p.jacobian.col(8).norm(), 0, 0, "Instantaneous trajectory independent of CR");
    // At a fixed raw clock tag, bias changes the reception geometry too.
    x[3] = 10;
    near(model(0, x).jacobian(0, 6), c - 10, .1, "Clock-dependent reception geometry");
    const fd::perturbations::AccelerationFunction srp = [](double, const auto &, const auto &) {
        fd::perturbations::AccelerationEvaluation a;
        a.acceleration[0] = 1e-7;
        return a;
    };
    const radio::TrajectoryFactory accelerated = [&](double epoch, const Eigen::VectorXd &s) {
        return std::make_unique<fd::dynamics::IntegratedStateProvider>(
            epoch, s.head<6>(), 101, [&, cr = s[8]](double t, const auto &r, const auto &v) {
                auto a = srp(t, r, v);
                a.acceleration *= cr;
                return a;
            });
    };
    const radio::OneWayMeasurement crModel(radio::OneWayLink(radio::libraryReceptionSolver()), tx,
                                           100, 60, 0, {0, f, 0}, f, {}, {}, 0, accelerated, steps);
    x[3] = 0;
    const auto crPrediction = crModel(0, x);
    near(crPrediction.jacobian(0, 8), .5 * 1e-7 * 10000, 1e-8, "Direct CR pseudorange partial");
    near(crPrediction.jacobian(1, 8), -f * 1e-7 * 70 / c, 1e-5,
         "Direct CR counted Doppler partial");
    fd::filters::CartesianPropagationConfig cfg;
    cfg.srpAcceleration = srp;
    cfg.clock = {0, 1e-24, 1e-30};
    fd::filters::CartesianPropagator propagate(srp, cfg);
    x[3] = 0;
    x[6] = 1e-6;
    x[7] = 2e-10;
    auto predicted = propagate(0, 100, x);
    near(predicted.state[0] - x[0], .5 * 1.2 * 1e-7 * 10000, 1e-8, "CR acceleration scaling");
    near(predicted.transition(0, 8), .5 * 1e-7 * 10000, 1e-10, "CR position STM");
    near(predicted.transition(3, 8), 1e-7 * 100, 1e-12, "CR velocity STM");
    near(predicted.state[6], 1e-6 + 2e-10 * 100, 1e-18, "Clock propagation");
    fd::filters::EKF filter(fd::filters::StateLayout::OrbitClockSrp);
    filter.setInitialState(x, steps.array().square().matrix().asDiagonal(), 0);
    filter.predictTo(100, propagate);
    const auto update =
        filter.update(100, model.value(100, filter.state()), Eigen::Matrix2d::Identity(),
                      [&](double t, const auto &s) { return model(t, s); });
    near(update.innovation.norm(), 0, 0, "Nine-state update");
}
} // namespace
int main() {
    try {
        physics();
        measurementAndNineStates();
        weakClockPartialAtDeepSpaceEpoch();
        std::cout << "One-way event, clocks, phase counts, group/phase corrections, Jacobians and "
                     "nine-state integration passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
