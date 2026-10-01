#pragma once

#include <Eigen/Core>

namespace fd::clocks {

// All times are elapsed seconds since calibration. b = clock - reference.
struct ClockState {
    double bias_s{0.0};
    double fractional_frequency{0.0};
};

struct ClockParameters {
    // SDE: db = y dt + sqrt(q_b) dW_b; dy = D dt + sqrt(q_y) dW_y.
    // Independent Wiener processes. These are diffusion intensities, not PSD
    // coefficients or per-step standard deviations; increments scale sqrt(dt).
    double frequency_drift_per_s{0.0};
    double q_bias_s{0.0};
    double q_frequency_per_s{0.0};
};

using ClockTransition = Eigen::Matrix2d;
using ClockCovariance = Eigen::Matrix2d;

} // namespace fd::clocks
