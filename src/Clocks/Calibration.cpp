#include "Clocks/Calibration.hpp"
#include "Validation.hpp"

#include <cmath>

namespace fd::clocks {

using namespace detail;

double clockRangeEnvelope(double t, const ClockParameters& p, const ClockBudget& b) {
    nonnegative(t);
    validate(p);
    validate(b.initial_state);
    nonnegative(b.sigma_bias_s);
    nonnegative(b.sigma_fractional_frequency);
    nonnegative(b.sigma_multiplier);
    const double deterministic = std::abs(b.initial_state.bias_s)
        + std::abs(b.initial_state.fractional_frequency) * t
        + (0.5 * std::abs(p.frequency_drift_per_s) * t) * t;
    const double sigma_y = b.sigma_fractional_frequency * t;
    const double variance = b.sigma_bias_s * b.sigma_bias_s + sigma_y * sigma_y
        + p.q_bias_s * t + ((p.q_frequency_per_s * t / 3.0) * t) * t;
    return checked(clockSpeedOfLightMPerS * (deterministic + b.sigma_multiplier * std::sqrt(variance)));
}

std::optional<double> clockCalibrationInterval(
    double threshold, double horizon, const ClockParameters& p, const ClockBudget& b) {
    positive(threshold);
    positive(horizon);
    if (clockRangeEnvelope(0.0, p, b) >= threshold) return 0.0;
    if (clockRangeEnvelope(horizon, p, b) < threshold) return std::nullopt;
    double lo = 0.0, hi = horizon;
    for (int i = 0; i < 100; ++i) {
        const double mid = lo + (hi - lo) / 2.0;
        if (mid == lo || mid == hi) break;
        if (clockRangeEnvelope(mid, p, b) < threshold) lo = mid;
        else hi = mid;
    }
    return hi;
}

} // namespace fd::clocks
