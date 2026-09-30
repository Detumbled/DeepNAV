#pragma once

#include <cstddef>
#include <vector>

namespace fd::clocks {

struct ClockSample {
    double time_s;
    double bias_s;
    double fractional_frequency;
    double sigma_bias_s{0.0};
    double sigma_fractional_frequency{0.0};
};

// One preallocated contiguous buffer; sigma fields are used only with covariance.
class ClockHistory {
public:
    explicit ClockHistory(std::size_t capacity, bool includeCovariance);
    ClockHistory(const ClockHistory&) = delete;
    ClockHistory& operator=(const ClockHistory&) = delete;
    ClockHistory(ClockHistory&&) noexcept = default;
    ClockHistory& operator=(ClockHistory&&) noexcept = default;

    [[nodiscard]] static std::size_t requiredCapacity(std::size_t steps, std::size_t recordEvery);
    [[nodiscard]] const std::vector<ClockSample>& samples() const noexcept { return samples_; }
    [[nodiscard]] bool includesCovariance() const noexcept { return includeCovariance_; }
    [[nodiscard]] std::size_t sampleLimit() const noexcept { return limit_; }

private:
    friend class Clocks;
    void record(const ClockSample& sample);
    void clear() noexcept { samples_.clear(); }
    std::vector<ClockSample> samples_;
    std::size_t limit_;
    bool includeCovariance_;
};

} // namespace fd::clocks
