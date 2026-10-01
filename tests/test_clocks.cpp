#include "Clocks/Allan.hpp"
#include "Clocks/Calibration.hpp"
#include "Clocks/ClockTruthSimulator.hpp"
#include "Clocks/DSAC.hpp"
#include "Clocks/LocalOscillator.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace fd::clocks;

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double relative = 1.0e-11, double absolute = 1.0e-30) {
    require(std::isfinite(actual) && std::abs(actual - expected) <= absolute + relative * std::abs(expected),
        "Numerical comparison failed");
}

template<class F> void rejects(F f) {
    try { f(); } catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("Invalid input was accepted");
}

void deterministicAndCovariance() {
    static_assert(std::is_abstract_v<Clocks> && std::has_virtual_destructor_v<Clocks>);
    const ClockParameters p{2e-15, 1e-24, 3e-26};
    const auto local = LocalOscillator::fromParameters(p);
    const auto dsac = DSAC::fromParameters(p);
    const Clocks& a = local;
    const Clocks& b = dsac;
    ClockState state{2e-9, 1e-11};
    ClockCovariance covariance;
    covariance << 4e-18, -1e-22, -1e-22, 9e-26;
    const auto initial_covariance = covariance;
    constexpr double dt = 7.0;
    for (int i = 0; i < 100; ++i) {
        const auto next = a.propagate(state, dt);
        const auto other = b.propagate(state, dt);
        near(next.bias_s, other.bias_s);
        near(next.fractional_frequency, other.fractional_frequency);
        require(a.transition(dt) == b.transition(dt), "Polymorphic F mismatch");
        require(a.processNoise(dt) == b.processNoise(dt), "Polymorphic Q mismatch");
        state = next;
        covariance = propagateClockCovariance(a, covariance, dt);
    }
    constexpr double t = 700.0;
    near(state.bias_s, 2e-9 + 1e-11 * t + 0.5 * p.frequency_drift_per_s * t * t);
    near(state.fractional_frequency, 1e-11 + p.frequency_drift_per_s * t);
    near(covariance(0, 0), initial_covariance(0, 0) + 2 * t * initial_covariance(0, 1)
        + t * t * initial_covariance(1, 1) + p.q_bias_s * t + p.q_frequency_per_s * t * t * t / 3);
    near(covariance(0, 1), initial_covariance(0, 1) + t * initial_covariance(1, 1)
        + p.q_frequency_per_s * t * t / 2);
    near(covariance(1, 1), initial_covariance(1, 1) + p.q_frequency_per_s * t);
    const auto q = a.processNoise(10.0);
    near(q(0, 0), 2e-23);
    near(q(0, 1), 1.5e-24);
    near(q(1, 1), 3e-25);
    std::cout << "PASS: deterministic propagation, polymorphism, F/Q and correlated initial covariance\n";
}

void allanFactories() {
    const ClockParameters p{2e-15, 1e-24, 3e-26};
    std::vector<AllanDatum> data;
    for (double tau : {1.0, 4.0, 16.0, 64.0, 256.0})
        data.push_back({tau, theoreticalAllanDeviation(tau, p)});
    AllanFitAssumptions assumptions{ClockNoiseAssumption::Both, true, 1.0, 256.0, 2e-15, 1e-10};
    const auto local = LocalOscillator::fromAllanData(data, assumptions);
    const auto dsac = DSAC::fromAllanData(data, assumptions);
    for (const Clocks* clock : {static_cast<const Clocks*>(&local), static_cast<const Clocks*>(&dsac)}) {
        near(clock->parameters().q_bias_s, p.q_bias_s);
        near(clock->parameters().q_frequency_per_s, p.q_frequency_per_s);
        near(clock->parameters().frequency_drift_per_s, p.frequency_drift_per_s);
    }
    const std::array<AllanDatum, 1> single{{{4.0, 2e-12}}};
    rejects([&] { (void)LocalOscillator::fromAllanData(single, assumptions); });
    rejects([&] { (void)DSAC::fromAllanData(single, assumptions); });
    assumptions.noise = ClockNoiseAssumption::WhiteFrequencyOnly;
    near(DSAC::fromAllanData(single, assumptions).parameters().q_bias_s, 16e-24);
    assumptions.noise = ClockNoiseAssumption::RandomWalkFrequencyOnly;
    near(LocalOscillator::fromAllanData(single, assumptions).parameters().q_frequency_per_s, 3e-24);
    assumptions.noise = ClockNoiseAssumption::Both;
    const std::array<AllanDatum, 2> repeated{{{4, 2e-12}, {4, 2e-12}}};
    rejects([&] { (void)fitClockAllanData(repeated, assumptions); });
    const std::array<AllanDatum, 2> close{{{4, 2e-12}, {4.0000000001, 2e-12}}};
    rejects([&] { (void)fitClockAllanData(close, assumptions); });
    const std::array<AllanDatum, 2> negative{{{1, 1e-12}, {4, 0.1e-12}}};
    rejects([&] { (void)fitClockAllanData(negative, assumptions); });
    const std::array<AllanDatum, 3> plateau{{{1, 1e-12}, {16, 1e-12}, {256, 1e-12}}};
    rejects([&] { (void)fitClockAllanData(plateau, assumptions); });
    assumptions.drift_removed = false;
    rejects([&] { (void)fitClockAllanData(data, assumptions); });
    assumptions.drift_removed = true;
    assumptions.valid_tau_min_s = 2;
    rejects([&] { (void)fitClockAllanData(data, assumptions); });
    assumptions.valid_tau_min_s = 1;
    for (auto& datum : data) datum.adev = 0;
    near(fitClockAllanData(data, assumptions).q_bias_s, 0);
    near(fitClockAllanData(data, assumptions).q_frequency_per_s, 0);
    std::cout << "PASS: Allan factories, dominant-noise conversion and incompatible/ambiguous fit rejection\n";
}

void referenceBaselines() {
    const auto clock = DSAC::dayMatchedWhiteFmBaseline();
    const auto& p = clock.parameters();
    require(clock.configurationName() == "DSAC_day_matched_white_FM_baseline", "Wrong baseline name");
    require(DSAC::fromParameters(p).configurationName() == "DSAC_custom",
        "Custom parameters must not implicitly claim a named baseline");
    near(p.frequency_drift_per_s, 3e-16 / 86400.0);
    near(p.q_bias_s, 7.776e-25);
    near(p.q_frequency_per_s, 0);
    near(theoreticalAllanDeviation(86400, p), 3e-15);
    const auto local = LocalOscillator::representativeUsoAgingOnly();
    near(local.parameters().frequency_drift_per_s, 1e-10 / 86400.0);
    near(local.parameters().q_bias_s, 0);
    near(local.parameters().q_frequency_per_s, 0);
    require(local.configurationName() == "USO_representative_aging_only", "Wrong USO reference name");

    const auto uso = LocalOscillator::representativeUsoWhiteFmWithAging();
    require(uso.configurationName() == "USO_aging_simplified_white_FM", "Wrong noisy USO name");
    near(uso.parameters().q_bias_s, 2.5e-25);
    near(uso.parameters().q_frequency_per_s, 0);
    near(uso.parameters().frequency_drift_per_s, 1e-10 / 86400.0);
    const auto uso_covariance = propagateClockCovariance(uso, ClockCovariance::Zero(), 20*86400);
    near(uso_covariance(0,0), 2.5e-25 * 20*86400);
    near(3 * clockSpeedOfLightMPerS * std::sqrt(uso_covariance(0,0)), 0.5911311305394723, 1e-10);
    near(clockSpeedOfLightMPerS * uso.propagate({}, 20*86400).bias_s, 518041.367424, 1e-10);
    for (const auto& row : {std::array<double,4>{2, 0.01554, 0.32968, 0.34522},
                           std::array<double,4>{10, 0.38853, 0.73719, 1.12572},
                           std::array<double,4>{20, 1.55412, 1.04254, 2.59666}}) {
        const double time = row[0]*86400;
        const double mean = clockSpeedOfLightMPerS * clock.propagate({}, time).bias_s;
        const double stochastic = 3*clockSpeedOfLightMPerS * std::sqrt(p.q_bias_s*time);
        require(std::abs(mean-row[1]) < 1e-5, "Incorrect DSAC drift contribution");
        require(std::abs(stochastic-row[2]) < 1e-5, "Incorrect DSAC stochastic contribution");
        require(std::abs(mean+stochastic-row[3]) < 1e-5, "Incorrect DSAC total budget");
    }
    constexpr double duration = 2 * 86400.0;
    for (double dt : {1.0, 60.0, 240.0, 86400.0}) {
        ClockCovariance covariance = ClockCovariance::Zero();
        for (int i = 0; i < static_cast<int>(duration / dt); ++i)
            covariance = propagateClockCovariance(clock, covariance, dt);
        near(covariance(0, 0), p.q_bias_s * duration, 1e-10);
        near(covariance(0, 1), 0);
        near(covariance(1, 1), 0);
        near(clock.processNoise(dt)(0, 0), p.q_bias_s * dt);
    }
    const double sigma = std::sqrt(p.q_bias_s * duration);
    near(sigma * 1e9, 0.367, 0.002);
    near(3 * clockSpeedOfLightMPerS * sigma, 0.330, 0.002);
    near(clockSpeedOfLightMPerS * clock.propagate({}, duration).bias_s, 0.0155, 0.003);

    const auto noise_free = DSAC::fromParameters({p.frequency_drift_per_s, 0, 0});
    ClockTruthSimulator noiseless(noise_free, 43);
    ClockState state{5e-9, 1e-12};
    for (int i = 0; i < 1000; ++i) state = noiseless.step(state, 7);
    constexpr double t = 7000;
    near(state.bias_s, 5e-9 + 1e-12*t + 0.5*p.frequency_drift_per_s*t*t);
    near(state.fractional_frequency, 1e-12 + p.frequency_drift_per_s*t);

    // Same seed gives the same standard normal at the first step. The white-FM
    // increment must scale as sqrt(dt), while y remains fully deterministic.
    ClockTruthSimulator reference(clock, 43);
    const double reference_noise = reference.step({}, 1).bias_s - clock.propagate({}, 1).bias_s;
    for (double dt : {4.0, 100.0, 86400.0}) {
        ClockTruthSimulator scaled(clock, 43);
        const auto draw = scaled.step({}, dt);
        const auto expected = clock.propagate({}, dt);
        near((draw.bias_s - expected.bias_s) / std::sqrt(dt), reference_noise, 1e-12);
        near(draw.fractional_frequency, expected.fractional_frequency);
    }
    ClockCovariance initial = ClockCovariance::Zero();
    initial(0, 0) = 4e-20;
    initial(1, 1) = 1e-30;
    const auto propagated = propagateClockCovariance(clock, initial, 100);
    near(propagated(0, 0), 4e-20 + 1e-26 + p.q_bias_s * 100);

    // A single long realization checks short tau only. It cannot validate
    // hardware day-scale ADEV from a two-day recording.
    ClockTruthSimulator truth(clock, 44);
    std::vector<double> bias(131073);
    ClockState current;
    for (std::size_t i = 1; i < bias.size(); ++i) {
        current = truth.step(current, 1);
        bias[i] = current.bias_s;
    }
    for (std::size_t m : {1, 4, 16, 64, 256})
        near(overlappingAllanDeviation(bias, 1, m) / std::sqrt(p.q_bias_s / m), 1, 0.20);
    std::cout << "PASS: named DSAC/USO baselines, two-day budget, step independence, sqrt(dt) noise and short-tau Allan\n";
}

void simulation() {
    // Synthetic coefficients only. Verify each noise source, including singular Q.
    for (const auto p : {ClockParameters{}, ClockParameters{0, 4e-24, 0},
                        ClockParameters{0, 0, 3e-26}, ClockParameters{0, 4e-24, 3e-26}}) {
        const auto clock = LocalOscillator::fromParameters(p);
        ClockTruthSimulator simulator(clock, 42);
        constexpr int count = 200000;
        Eigen::Vector2d sum = Eigen::Vector2d::Zero();
        Eigen::Matrix2d sum_products = Eigen::Matrix2d::Zero();
        for (int i = 0; i < count; ++i) {
            const auto s = simulator.step({}, 10);
            const Eigen::Vector2d v{s.bias_s, s.fractional_frequency};
            sum += v;
            sum_products += v * v.transpose();
        }
        const Eigen::Matrix2d empirical = (sum_products - sum * sum.transpose() / count) / (count - 1);
        const auto q = clock.processNoise(10);
        for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c) near(empirical(r, c), q(r, c), 0.03);
        ClockTruthSimulator other(clock, 42);
        simulator.reseed(42);
        for (int i = 0; i < 10; ++i) {
            const auto x = simulator.step({}, 1);
            const auto y = other.step({}, 1);
            require(x.bias_s == y.bias_s && x.fractional_frequency == y.fractional_frequency,
                "Reseeding does not reproduce samples");
        }
    }
    std::cout << "PASS: empirical covariance including cross terms, zero/single noise and reseeding\n";
}

void allanStatistics() {
    int seed = 100;
    for (const auto p : {ClockParameters{0, 1e-24, 0}, ClockParameters{0, 0, 3e-26},
                        ClockParameters{0, 1e-24, 3e-26}}) {
        const auto clock = DSAC::fromParameters(p);
        ClockTruthSimulator simulator(clock, seed++);
        std::vector<double> bias(262145);
        ClockState state;
        for (std::size_t i = 1; i < bias.size(); ++i) {
            state = simulator.step(state, 1);
            bias[i] = state.bias_s;
        }
        std::cout << "PASS: simulated/theoretical ADEV";
        for (std::size_t m : {1, 4, 16, 64, 256}) {
            const double ratio = overlappingAllanDeviation(bias, 1, m) / theoreticalAllanDeviation(m, p);
            near(ratio, 1, 0.20);
            std::cout << ' ' << m << "s:" << ratio;
        }
        std::cout << '\n';
    }
    std::vector<double> linear(4097), quadratic(linear.size());
    for (std::size_t i = 0; i < linear.size(); ++i) {
        linear[i] = 1e-8 + 2e-11 * i;
        quadratic[i] = 0.5 * 2e-15 * i * i;
    }
    require(overlappingAllanDeviation(linear, 1, 64) < 1e-22, "ADEV sees constant bias/frequency");
    near(overlappingAllanDeviation(quadratic, 1, 64), 2e-15 * 64 / std::sqrt(2.0), 1e-9);
    near(theoreticalAllanDeviation(64, {2e-15, 0, 0}, true), 2e-15 * 64 / std::sqrt(2.0));
    near(theoreticalAllanDeviation(64, {2e-15, 0, 0}), 0);
    std::cout << "PASS: Allan invariance and unremoved drift\n";
}

void calibration() {
    const double c = clockSpeedOfLightMPerS, k = 3;
    constexpr double horizon = 1e9;
    ClockBudget budget;
    budget.initial_state.fractional_frequency = -3e-15;
    near(*clockCalibrationInterval(1, horizon, {}, budget), 1 / (c * 3e-15));
    budget = {};
    near(*clockCalibrationInterval(1, horizon, {2e-15, 0, 0}, budget), std::sqrt(2 / (c * 2e-15)));
    near(*clockCalibrationInterval(1, horizon, {0, 1e-24, 0}, budget), std::pow(1 / (c * k), 2) / 1e-24);
    near(*clockCalibrationInterval(1, horizon, {0, 0, 3e-26}, budget), std::cbrt(3 * std::pow(1 / (c * k), 2) / 3e-26));
    require(!clockCalibrationInterval(1, 86400, {}), "Ideal clock crossed budget");
    budget.initial_state.bias_s = 1e-8;
    near(*clockCalibrationInterval(1, 86400, {}, budget), 0);
    budget.initial_state.bias_s = 1 / c;
    near(*clockCalibrationInterval(1, 86400, {}, budget), 0);
    budget = {};
    budget.sigma_bias_s = 2e-10;
    budget.sigma_fractional_frequency = 1e-15;
    near(clockRangeEnvelope(10, {0, 0, 0}, budget), c * k * std::sqrt(4e-20 + 1e-28));
    // Deterministic signs cannot cancel and extend a calibration interval.
    budget = {};
    budget.initial_state = {1e-9, -1e-12};
    near(clockRangeEnvelope(1000, {0, 0, 0}, budget), c * 2e-9);
    std::cout << "PASS: four analytic calibration limits, initial contact, horizon and conservative envelope\n";
}

void validation() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    rejects([&] { (void)LocalOscillator::fromParameters({nan, 0, 0}); });
    rejects([&] { (void)DSAC::fromParameters({0, -1, 0}); });
    rejects([&] { (void)DSAC::fromParameters({0, 0, inf}); });
    const auto clock = DSAC::fromParameters({});
    for (double dt : {0.0, -1.0, nan, inf}) {
        rejects([&] { (void)clock.propagate({}, dt); });
        rejects([&] { (void)clock.transition(dt); });
        rejects([&] { (void)clock.processNoise(dt); });
    }
    rejects([&] { (void)clock.propagate({nan, 0}, 1); });
    const std::array<double, 3> bias{0, 1, 2};
    rejects([&] { (void)overlappingAllanDeviation(bias, 1, 0); });
    rejects([&] { (void)overlappingAllanDeviation(bias, 1, 2); });
    rejects([&] { (void)overlappingAllanDeviation(bias, 1, std::numeric_limits<std::size_t>::max()); });
    near(overlappingAllanDeviation(bias, 1, 1), 0);
    const std::array<double, 3> bad_bias{nan, 0, 0};
    rejects([&] { (void)overlappingAllanDeviation(bad_bias, 1, 1); });
    rejects([&] { (void)theoreticalAllanDeviation(0, {}); });
    rejects([&] { (void)clockRangeEnvelope(-1, {}); });
    rejects([&] { (void)clockCalibrationInterval(0, 1, {}); });
    rejects([&] { (void)clockCalibrationInterval(1, inf, {}); });
    ClockBudget bad_budget;
    bad_budget.sigma_bias_s = -1;
    rejects([&] { (void)clockRangeEnvelope(1, {}, bad_budget); });
    ClockCovariance bad_covariance;
    bad_covariance << 1, 2, 2, 1;
    rejects([&] { (void)propagateClockCovariance(clock, bad_covariance, 1); });
    std::cout << "PASS: invalid inputs rejected\n";
}

} // namespace

int main() {
    try {
        deterministicAndCovariance();
        allanFactories();
        referenceBaselines();
        simulation();
        allanStatistics();
        calibration();
        validation();
        std::cout << "All clock tests passed (synthetic tests and documented effective baselines).\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
