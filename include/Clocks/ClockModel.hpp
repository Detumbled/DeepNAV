#pragma once

#include "Clocks/Types.hpp"

namespace fd::clocks {

// Shared two-state SDE: white frequency and random walk frequency noise only.
// Diffusion intensities are not one-sided PSD coefficients.
class ClockModel {
public:
    explicit ClockModel(ClockParameters parameters);
    [[nodiscard]] ClockState propagate(const ClockState& state, double dt) const;
    [[nodiscard]] ClockTransition transition(double dt) const;
    [[nodiscard]] ClockCovariance processNoise(double dt) const;
    [[nodiscard]] const ClockParameters& parameters() const noexcept { return parameters_; }
private:
    ClockParameters parameters_;
};

} // namespace fd::clocks
