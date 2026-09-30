#include "Clocks/LocalOscillator.hpp"

namespace fd::clocks {

LocalOscillator LocalOscillator::fromParameters(ClockParameters p) { return LocalOscillator(p); }
LocalOscillator LocalOscillator::fromAllanData(std::span<const AllanDatum> d, const AllanFitAssumptions& a) {
    return fromParameters(fitClockAllanData(d, a));
}
ClockState LocalOscillator::propagate(const ClockState& s, double dt) const { return model_.propagate(s, dt); }
ClockTransition LocalOscillator::transition(double dt) const { return model_.transition(dt); }
ClockCovariance LocalOscillator::processNoise(double dt) const { return model_.processNoise(dt); }
const ClockParameters& LocalOscillator::parameters() const noexcept { return model_.parameters(); }

} // namespace fd::clocks
