#pragma once

#include <Eigen/Core>

namespace fd::clocks {

// All times are elapsed seconds since calibration. b = clock - reference.
struct ClockState {
    double bias_s{0.0};
    double fractional_frequency{0.0};
};

struct ClockParameters {
    double frequency_drift_per_s{0.0};
    double q_bias_s{0.0};
    double q_frequency_per_s{0.0};
};

using ClockTransition = Eigen::Matrix2d;
using ClockCovariance = Eigen::Matrix2d;

} // namespace fd::clocks
