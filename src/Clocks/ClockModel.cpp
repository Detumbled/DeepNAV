#include "Clocks/ClockModel.hpp"
#include "Validation.hpp"

namespace fd::clocks {

using namespace detail;

ClockModel::ClockModel(ClockParameters parameters) : parameters_(parameters) {
    validate(parameters_);
}

ClockState ClockModel::propagate(const ClockState& state, double dt) const {
    positive(dt);
    validate(state);
    const double drift_step = checked(parameters_.frequency_drift_per_s * dt);
    return {checked(state.bias_s + state.fractional_frequency * dt + 0.5 * drift_step * dt),
            checked(state.fractional_frequency + drift_step)};
}

ClockTransition ClockModel::transition(double dt) const {
    positive(dt);
    ClockTransition f;
    f << 1.0, dt, 0.0, 1.0;
    return f;
}

ClockCovariance ClockModel::processNoise(double dt) const {
    positive(dt);
    const double yy = checked(parameters_.q_frequency_per_s * dt);
    const double by = checked(0.5 * yy * dt);
    const double bb = checked(parameters_.q_bias_s * dt + (yy * dt / 3.0) * dt);
    ClockCovariance q;
    q << bb, by, by, yy;
    return q;
}

} // namespace fd::clocks
