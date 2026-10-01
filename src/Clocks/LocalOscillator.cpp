#include "Clocks/LocalOscillator.hpp"

namespace fd::clocks {

LocalOscillator LocalOscillator::representativeUsoAgingOnly() {
    // doi:10.1029/2025RS008244, Table 2. Aging is fractional frequency/day.
    // Two tabulated ADEV points are not used to invent q_b or q_y.
    return LocalOscillator({1.0e-10 / 86400.0, 0.0, 0.0}, "USO_representative_aging_only");
}

LocalOscillator LocalOscillator::representativeUsoWhiteFmWithAging() {
    return LocalOscillator({1e-10 / 86400.0, 2.5e-25, 0.0},
                           "USO_aging_simplified_white_FM");
}

LocalOscillator LocalOscillator::fromParameters(ClockParameters p) { return LocalOscillator(p); }
LocalOscillator LocalOscillator::fromAllanData(std::span<const AllanDatum> d, const AllanFitAssumptions& a) {
    return fromParameters(fitClockAllanData(d, a));
}
ClockState LocalOscillator::propagate(const ClockState& s, double dt) const { return model_.propagate(s, dt); }
ClockTransition LocalOscillator::transition(double dt) const { return model_.transition(dt); }
ClockCovariance LocalOscillator::processNoise(double dt) const { return model_.processNoise(dt); }
const ClockParameters& LocalOscillator::parameters() const noexcept { return model_.parameters(); }

} // namespace fd::clocks
