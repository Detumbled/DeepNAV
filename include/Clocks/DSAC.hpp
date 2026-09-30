#pragma once

#include "Clocks/Allan.hpp"
#include "Clocks/ClockModel.hpp"
#include "Clocks/Clocks.hpp"

namespace fd::clocks {

// Both factories produce custom configurations, not certified DSAC presets.
class DSAC final : public Clocks {
public:
    [[nodiscard]] static DSAC fromParameters(ClockParameters parameters);
    [[nodiscard]] static DSAC fromAllanData(
        std::span<const AllanDatum> data, const AllanFitAssumptions& assumptions);
    [[nodiscard]] ClockState propagate(const ClockState& state, double dt) const override;
    [[nodiscard]] ClockTransition transition(double dt) const override;
    [[nodiscard]] ClockCovariance processNoise(double dt) const override;
    [[nodiscard]] const ClockParameters& parameters() const noexcept override;
private:
    explicit DSAC(ClockParameters parameters) : model_(parameters) {}
    ClockModel model_;
};

} // namespace fd::clocks
