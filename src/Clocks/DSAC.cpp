#include "Clocks/DSAC.hpp"

namespace fd::clocks {

DSAC DSAC::shortTermWhiteFmBaseline() {
    // Burt et al. (2021), doi:10.1038/s41586-021-03571-7, abstract:
    // ground-test short-term ADEV 1.5e-13/sqrt(tau/1 s).
    // Known flight drift is separate; no long-term stochastic fit is implied.
    constexpr double shortTermAdev = 1.5e-13;
    return DSAC({3e-16 / 86400.0, shortTermAdev * shortTermAdev, 0.0},
                "DSAC_inspired_short_term_white_FM");
}

DSAC DSAC::dayMatchedWhiteFmBaseline() {
    // Tjoelker (2021), slide 21: 3e-15 one-day ADEV.
    // Burt et al. (2021), doi:10.1038/s41586-021-03571-7: 3e-16/day drift.
    // The q_b conversion is an effective white-FM approximation, not a fit.
    constexpr double secondsPerDay = 86400.0;
    constexpr double dayAdev = 3.0e-15;
    return DSAC({3.0e-16 / secondsPerDay, dayAdev * dayAdev * secondsPerDay, 0.0},
                "DSAC_day_matched_white_FM_baseline");
}

DSAC DSAC::fromParameters(ClockParameters p) { return DSAC(p); }
DSAC DSAC::fromAllanData(std::span<const AllanDatum> d, const AllanFitAssumptions& a) {
    return fromParameters(fitClockAllanData(d, a));
}
ClockState DSAC::propagate(const ClockState& s, double dt) const { return model_.propagate(s, dt); }
ClockTransition DSAC::transition(double dt) const { return model_.transition(dt); }
ClockCovariance DSAC::processNoise(double dt) const { return model_.processNoise(dt); }
const ClockParameters& DSAC::parameters() const noexcept { return model_.parameters(); }

} // namespace fd::clocks
