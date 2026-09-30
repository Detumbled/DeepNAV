#pragma once

#include "Clocks/Types.hpp"
#include "Clocks/ClockHistory.hpp"

#include <memory>

namespace fd::clocks {

class Clocks {
public:
    Clocks() = default;
    Clocks(const Clocks&) = delete;
    Clocks& operator=(const Clocks&) = delete;
    Clocks(Clocks&&) noexcept = default;
    Clocks& operator=(Clocks&&) noexcept = default;
    virtual ~Clocks() = default;
    [[nodiscard]] virtual ClockState propagate(const ClockState& state, double dt) const = 0;
    [[nodiscard]] virtual ClockTransition transition(double dt) const = 0;
    [[nodiscard]] virtual ClockCovariance processNoise(double dt) const = 0;
    [[nodiscard]] virtual const ClockParameters& parameters() const noexcept = 0;

    void enableHistory(std::size_t capacity, bool includeCovariance);
    void disableHistory() noexcept;
    void clearHistory() noexcept;
    [[nodiscard]] bool historyEnabled() const noexcept { return history_ != nullptr; }
    [[nodiscard]] const ClockHistory* history() const noexcept { return history_.get(); }
    // Explicit recording only: propagate() and simulator.step() remain stateless.
    void recordSample(const ClockSample& sample);

private:
    std::unique_ptr<ClockHistory> history_;
};

[[nodiscard]] ClockCovariance propagateClockCovariance(
    const Clocks& clock, const ClockCovariance& covariance, double dt);

} // namespace fd::clocks
