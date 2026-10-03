#pragma once

#include "Clocks/Allan.hpp"
#include "Clocks/ClockModel.hpp"
#include "Clocks/Clocks.hpp"

#include <string_view>

namespace fd::clocks {

// Custom configurations and a documented effective baseline; no hardware fit.
class DSAC final : public Clocks {
public:
    // DSAC-inspired short-term ground-test law: ADEV = 1.5e-13/sqrt(tau/1 s).
    // q_b=2.25e-26 s, q_y=0; deterministic drift remains 3e-16/day.
    // Long-term white-FM continuation is illustrative, not a hardware fit.
    [[nodiscard]] static DSAC shortTermWhiteFmBaseline();
    // Legacy one-day approximation, retained for explicit comparisons only.
    // White FM matched to stochastic ADEV 3e-15 at 86400 s, with known drift
    // 3e-16/day. Initial bias/frequency and P0 remain caller-owned (default zero).
    // Does not reproduce DSAC's short-term behavior or long-term floor.
    [[nodiscard]] static DSAC dayMatchedWhiteFmBaseline();
    [[nodiscard]] static DSAC fromParameters(ClockParameters parameters);
    [[nodiscard]] static DSAC fromAllanData(
        std::span<const AllanDatum> data, const AllanFitAssumptions& assumptions);
    [[nodiscard]] ClockState propagate(const ClockState& state, double dt) const override;
    [[nodiscard]] ClockTransition transition(double dt) const override;
    [[nodiscard]] ClockCovariance processNoise(double dt) const override;
    [[nodiscard]] const ClockParameters& parameters() const noexcept override;
    [[nodiscard]] std::string_view configurationName() const noexcept { return configurationName_; }
private:
    explicit DSAC(ClockParameters parameters, std::string_view name = "DSAC_custom")
        : model_(parameters), configurationName_(name) {}
    ClockModel model_;
    std::string_view configurationName_;
};

} // namespace fd::clocks
