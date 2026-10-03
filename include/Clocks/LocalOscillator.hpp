#pragma once

#include "Clocks/Allan.hpp"
#include "Clocks/ClockModel.hpp"
#include "Clocks/Clocks.hpp"

#include <string_view>

namespace fd::clocks {

class LocalOscillator final : public Clocks {
public:
    // Ely et al. (2025), Table 2: representative USO aging 1e-10/day.
    // Deterministic comparison only: no random noise fit or hardware prediction.
    [[nodiscard]] static LocalOscillator representativeUsoAgingOnly();
    // Representative USO-like illustrative white FM: ADEV 5e-13 at 1 s,
    // plus uncompensated aging 1e-10/day. Not a particular device noise fit.
    [[nodiscard]] static LocalOscillator representativeUsoWhiteFmWithAging();
    [[nodiscard]] static LocalOscillator fromParameters(ClockParameters parameters);
    [[nodiscard]] static LocalOscillator fromAllanData(
        std::span<const AllanDatum> data, const AllanFitAssumptions& assumptions);
    [[nodiscard]] ClockState propagate(const ClockState& state, double dt) const override;
    [[nodiscard]] ClockTransition transition(double dt) const override;
    [[nodiscard]] ClockCovariance processNoise(double dt) const override;
    [[nodiscard]] const ClockParameters& parameters() const noexcept override;
    [[nodiscard]] std::string_view configurationName() const noexcept { return configurationName_; }
private:
    explicit LocalOscillator(ClockParameters parameters, std::string_view name = "LocalOscillator_custom")
        : model_(parameters), configurationName_(name) {}
    ClockModel model_;
    std::string_view configurationName_;
};

} // namespace fd::clocks
