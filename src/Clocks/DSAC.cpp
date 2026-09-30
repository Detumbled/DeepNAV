#include "Clocks/DSAC.hpp"

namespace fd::clocks {

DSAC DSAC::fromParameters(ClockParameters p) { return DSAC(p); }
DSAC DSAC::fromAllanData(std::span<const AllanDatum> d, const AllanFitAssumptions& a) {
    return fromParameters(fitClockAllanData(d, a));
}
ClockState DSAC::propagate(const ClockState& s, double dt) const { return model_.propagate(s, dt); }
ClockTransition DSAC::transition(double dt) const { return model_.transition(dt); }
ClockCovariance DSAC::processNoise(double dt) const { return model_.processNoise(dt); }
const ClockParameters& DSAC::parameters() const noexcept { return model_.parameters(); }

} // namespace fd::clocks
