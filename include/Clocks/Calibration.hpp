#pragma once

#include "Clocks/Types.hpp"

#include <optional>

namespace fd::clocks {

inline constexpr double clockSpeedOfLightMPerS = 299792458.0;

struct ClockBudget {
    ClockState initial_state{};
    double sigma_bias_s{0.0};
    double sigma_fractional_frequency{0.0};
    double sigma_multiplier{3.0};
};

// Monotonic conservative envelope, assuming independent initial uncertainties.
// This pointwise Gaussian budget does not bound an entire sample path.
[[nodiscard]] double clockRangeEnvelope(
    double elapsed_s, const ClockParameters& parameters, const ClockBudget& budget = {});
// nullopt means > horizon_s, never infinite autonomy. Initial contact returns 0.
[[nodiscard]] std::optional<double> clockCalibrationInterval(
    double threshold_m, double horizon_s, const ClockParameters& parameters,
    const ClockBudget& budget = {});

} // namespace fd::clocks
