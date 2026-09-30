#pragma once

#include "Clocks/Allan.hpp"
#include "Clocks/ClockModel.hpp"
#include "Clocks/Clocks.hpp"

namespace fd::clocks {

class LocalOscillator final : public Clocks {
public:
    [[nodiscard]] static LocalOscillator fromParameters(ClockParameters parameters);
    [[nodiscard]] static LocalOscillator fromAllanData(
        std::span<const AllanDatum> data, const AllanFitAssumptions& assumptions);
    [[nodiscard]] ClockState propagate(const ClockState& state, double dt) const override;
    [[nodiscard]] ClockTransition transition(double dt) const override;
    [[nodiscard]] ClockCovariance processNoise(double dt) const override;
    [[nodiscard]] const ClockParameters& parameters() const noexcept override;
private:
    explicit LocalOscillator(ClockParameters parameters) : model_(parameters) {}
    ClockModel model_;
};

} // namespace fd::clocks
