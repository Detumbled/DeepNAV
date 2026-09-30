#include "Clocks/Clocks.hpp"
#include "Validation.hpp"

#include <cmath>
#include <stdexcept>

namespace fd::clocks {

void Clocks::enableHistory(std::size_t capacity, bool includeCovariance) {
    history_ = std::make_unique<ClockHistory>(capacity, includeCovariance);
}

void Clocks::disableHistory() noexcept { history_.reset(); }

void Clocks::clearHistory() noexcept {
    if (history_) history_->clear();
}

void Clocks::recordSample(const ClockSample& sample) {
    if (history_) history_->record(sample);
}

ClockCovariance propagateClockCovariance(const Clocks& clock, const ClockCovariance& p, double dt) {
    detail::validate(p);
    const auto f = clock.transition(dt);
    const ClockCovariance next = f * p * f.transpose() + clock.processNoise(dt);
    if (!next.allFinite()) throw std::overflow_error("Clock covariance overflowed.");
    // Preserve exact symmetry across repeated matrix products.
    return 0.5 * next + 0.5 * next.transpose();
}

} // namespace fd::clocks
