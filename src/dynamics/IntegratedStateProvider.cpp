#include "dynamics/IntegratedStateProvider.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::dynamics {
IntegratedStateProvider::IntegratedStateProvider(
    double initialTdb, const Eigen::Matrix<double, 6, 1> &initial, double endTdb,
    const perturbations::AccelerationFunction &acceleration, od::RKF45Integrator::Options options,
    std::string frame, std::string origin, const StateProvider *movingOrigin)
    : initialTdb_(initialTdb), frame_(std::move(frame)), origin_(std::move(origin)),
      movingOrigin_(movingOrigin) {
    if (!acceleration || frame_.empty() || origin_.empty() || !std::isfinite(initialTdb) ||
        !std::isfinite(endTdb) || endTdb <= initialTdb || !initial.allFinite())
        throw std::invalid_argument("Invalid numerical trajectory configuration.");
    if (movingOrigin_ &&
        (movingOrigin_->referenceFrame() != frame_ || movingOrigin_->referenceOrigin() != origin_))
        throw std::invalid_argument("Moving origin must share output frame and fixed origin.");
    const auto derivative = [&](double elapsed, od::RKF45Integrator::ConstStateRef x,
                                od::RKF45Integrator::StateRef rate) {
        const auto force = acceleration(initialTdb + elapsed, x.head<3>(), x.tail<3>());
        rate.head<3>() = x.tail<3>();
        rate.tail<3>() = force.acceleration;
    };
    // Relative integration times preserve small steps at large absolute epochs.
    const auto integrated =
        od::RKF45Integrator(options).integrate(0, initial, endTdb - initialTdb, derivative);
    for (const auto &node : integrated.history)
        ephemeris_.addNode(node.tdb, node.state, node.derivative);
}
CartesianState IntegratedStateProvider::relativeStateAt(double epoch) const {
    const auto x = ephemeris_.interpolate(epoch - initialTdb_);
    return {x.head<3>(), x.tail<3>()};
}
CartesianState IntegratedStateProvider::stateAt(opnav::TdbEpoch epoch) const {
    auto state = relativeStateAt(epoch.secondsPastJ2000);
    if (movingOrigin_) {
        const auto origin = movingOrigin_->stateAt(epoch);
        if (!origin.allFinite())
            throw std::runtime_error("Non-finite moving-origin state.");
        state.positionKm += origin.positionKm;
        state.velocityKmPerSec += origin.velocityKmPerSec;
    }
    return state;
}
} // namespace fd::dynamics
