#pragma once

#include "Clocks/Types.hpp"

#include <cstddef>
#include <span>

namespace fd::clocks {

struct AllanDatum {
    double tau_s;
    double adev;
};

enum class ClockNoiseAssumption { WhiteFrequencyOnly, RandomWalkFrequencyOnly, Both };

struct AllanFitAssumptions {
    ClockNoiseAssumption noise{ClockNoiseAssumption::Both};
    // Caller must explicitly confirm zero drift or removal before fitting.
    bool drift_removed{false};
    double valid_tau_min_s{0.0};
    double valid_tau_max_s{0.0};
    double frequency_drift_per_s{0.0}; // Known D, supplied separately from ADEV.
    double max_relative_variance_residual{0.05};
};

// Fits Allan variance, rejecting ill-conditioned, negative or incompatible fits.
[[nodiscard]] ClockParameters fitClockAllanData(
    std::span<const AllanDatum> data, const AllanFitAssumptions& assumptions);

[[nodiscard]] double overlappingAllanDeviation(
    std::span<const double> bias_s, double dt, std::size_t m);
[[nodiscard]] double theoreticalAllanDeviation(
    double tau_s, const ClockParameters& parameters, bool include_drift = false);

} // namespace fd::clocks
