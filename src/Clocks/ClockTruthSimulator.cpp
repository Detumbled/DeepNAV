#include "Clocks/ClockTruthSimulator.hpp"
#include "Validation.hpp"

#include <cmath>

namespace fd::clocks {

using namespace detail;

ClockTruthSimulator::ClockTruthSimulator(const Clocks& clock, std::uint64_t seed)
    : clock_(clock), generator_(seed) {}

ClockTruthSimulator::ClockTruthSimulator(Clocks& clock, std::uint64_t seed)
    : clock_(clock), recordingClock_(&clock), generator_(seed) {}

void ClockTruthSimulator::prepareRun(std::size_t steps, std::size_t recordEvery,
    const ClockState& initial, const std::optional<ClockCovariance>& covariance) const {
    validate(initial);
    const auto required = ClockHistory::requiredCapacity(steps, recordEvery);
    if (covariance) validate(*covariance);
    if (const auto* history = clock_.history()) {
        if (!recordingClock_)
            throw std::invalid_argument("Recording a run requires a mutable clock reference.");
        if (!history->samples().empty())
            throw std::invalid_argument("Clear clock history before starting a new run.");
        if (history->sampleLimit() < required)
            throw std::length_error("Clock history capacity is too small for this run and recordEvery.");
        if (history->includesCovariance() && !covariance)
            throw std::invalid_argument("Covariance history requires an explicit initial covariance.");
    }
}

void ClockTruthSimulator::record(double elapsed, const ClockState& state,
    const std::optional<ClockCovariance>& covariance) {
    if (!clock_.historyEnabled()) return;
    ClockSample sample{elapsed, state.bias_s, state.fractional_frequency};
    if (clock_.history()->includesCovariance()) {
        sample.sigma_bias_s = std::sqrt((*covariance)(0, 0));
        sample.sigma_fractional_frequency = std::sqrt((*covariance)(1, 1));
    }
    recordingClock_->recordSample(sample);
}

ClockState ClockTruthSimulator::run(std::size_t steps, double dt, ClockState state,
    std::size_t recordEvery, std::optional<ClockCovariance> covariance) {
    positive(dt);
    prepareRun(steps, recordEvery, state, covariance);
    double elapsed = 0.0;
    record(elapsed, state, covariance);
    for (std::size_t i = 0; i < steps; ++i) {
        const double next_time = checked(elapsed + dt);
        if (next_time <= elapsed) throw std::overflow_error("Clock time step cannot advance elapsed time.");
        state = step(state, dt);
        if (covariance) *covariance = propagateClockCovariance(clock_, *covariance, dt);
        elapsed = next_time;
        if ((i + 1) % recordEvery == 0 || i + 1 == steps) record(elapsed, state, covariance);
    }
    return state;
}

ClockState ClockTruthSimulator::run(std::span<const double> timeSteps, ClockState state,
    std::size_t recordEvery, std::optional<ClockCovariance> covariance) {
    for (double dt : timeSteps) positive(dt);
    prepareRun(timeSteps.size(), recordEvery, state, covariance);
    double elapsed = 0.0;
    record(elapsed, state, covariance);
    for (std::size_t i = 0; i < timeSteps.size(); ++i) {
        const double next_time = checked(elapsed + timeSteps[i]);
        if (next_time <= elapsed) throw std::overflow_error("Clock time step cannot advance elapsed time.");
        state = step(state, timeSteps[i]);
        if (covariance) *covariance = propagateClockCovariance(clock_, *covariance, timeSteps[i]);
        elapsed = next_time;
        if ((i + 1) % recordEvery == 0 || i + 1 == timeSteps.size()) record(elapsed, state, covariance);
    }
    return state;
}

void ClockTruthSimulator::reseed(std::uint64_t seed) {
    generator_.seed(seed);
    normal_.reset();
}

ClockState ClockTruthSimulator::step(const ClockState& state, double dt) {
    auto next = clock_.propagate(state, dt);
    const auto q = clock_.processNoise(dt); // Also rejects overflow before drawing.
    const double z0 = normal_(generator_);
    const double z1 = normal_(generator_);
    const double z2 = normal_(generator_);
    const double white = checked(clock_.parameters().q_bias_s * dt);
    const double integrated_walk = checked((q(1, 1) * dt) * dt);
    next.bias_s = checked(next.bias_s + std::sqrt(white) * z0
        + std::sqrt(integrated_walk) * (z1 / 2.0 + z2 / std::sqrt(12.0)));
    next.fractional_frequency = checked(next.fractional_frequency + std::sqrt(q(1, 1)) * z1);
    return next;
}

} // namespace fd::clocks
