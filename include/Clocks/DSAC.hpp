#pragma once

#include "Clocks/Allan.hpp"
#include "Clocks/ClockModel.hpp"
#include "Clocks/Clocks.hpp"

#include <string_view>

namespace fd::clocks {

// Custom configurations and a documented effective baseline; no hardware fit.
class DSAC final : public Clocks {
public:
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
