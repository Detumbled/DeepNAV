#include "Clocks/ClockHistory.hpp"
#include "Validation.hpp"

#include <limits>
#include <stdexcept>

namespace fd::clocks {

ClockHistory::ClockHistory(std::size_t capacity, bool includeCovariance)
    : limit_(capacity), includeCovariance_(includeCovariance) {
    samples_.reserve(capacity);
}

std::size_t ClockHistory::requiredCapacity(std::size_t steps, std::size_t recordEvery) {
    if (recordEvery == 0) throw std::invalid_argument("recordEvery must be >= 1.");
    const std::size_t extra = 1 + (steps % recordEvery != 0);
    const std::size_t selected = steps / recordEvery;
    if (selected > std::numeric_limits<std::size_t>::max() - extra)
        throw std::overflow_error("Clock history capacity overflow.");
    return selected + extra;
}

void ClockHistory::record(const ClockSample& sample) {
    detail::nonnegative(sample.time_s);
    detail::finite(sample.bias_s);
    detail::finite(sample.fractional_frequency);
    if (includeCovariance_) {
        detail::nonnegative(sample.sigma_bias_s);
        detail::nonnegative(sample.sigma_fractional_frequency);
    }
    if (!samples_.empty() && sample.time_s <= samples_.back().time_s)
        throw std::invalid_argument("Clock history timestamps must be strictly increasing.");
    if (samples_.size() >= limit_)
        throw std::length_error("Clock history sample limit reached; reserve a larger history before running.");
    samples_.push_back(sample);
}

} // namespace fd::clocks
