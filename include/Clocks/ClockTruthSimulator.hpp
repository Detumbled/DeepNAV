#pragma once

#include "Clocks/Clocks.hpp"

#include <cstdint>
#include <optional>
#include <random>
#include <span>

namespace fd::clocks {

class ClockTruthSimulator {
public:
    // clock must outlive the simulator. Truth and covariance are kept separate.
    explicit ClockTruthSimulator(const Clocks& clock, std::uint64_t seed = 2026);
    explicit ClockTruthSimulator(Clocks& clock, std::uint64_t seed = 2026);
    ClockTruthSimulator(const Clocks&&, std::uint64_t = 2026) = delete;
    [[nodiscard]] ClockState step(const ClockState& state, double dt);
    // A fresh run starts at elapsed time zero and requires an empty history.
    // Enable history on a mutable clock first; covariance history requires P0.
    [[nodiscard]] ClockState run(std::size_t steps, double dt, ClockState initial = {},
        std::size_t recordEvery = 1, std::optional<ClockCovariance> initialCovariance = std::nullopt);
    // Variable step sizes use actual accumulated elapsed time and an explicit
    // history capacity supplied by the caller via enableHistory().
    [[nodiscard]] ClockState run(std::span<const double> timeSteps, ClockState initial = {},
        std::size_t recordEvery = 1, std::optional<ClockCovariance> initialCovariance = std::nullopt);
    void reseed(std::uint64_t seed);
private:
    void prepareRun(std::size_t steps, std::size_t recordEvery, const ClockState& initial,
        const std::optional<ClockCovariance>& initialCovariance) const;
    void record(double elapsed, const ClockState& state,
        const std::optional<ClockCovariance>& covariance);
    const Clocks& clock_;
    Clocks* recordingClock_{nullptr};
    std::mt19937_64 generator_;
    std::normal_distribution<double> normal_{0.0, 1.0};
};

} // namespace fd::clocks
