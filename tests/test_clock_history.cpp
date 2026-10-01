#include "Clocks/ClockTruthSimulator.hpp"
#include "Clocks/DSAC.hpp"
#include "Clocks/LocalOscillator.hpp"
#include "utils/ClockCsvWriter.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

using namespace fd::clocks;
using fd::utils::writeClockCsv;
using fd::utils::ClockCsvWriteMode;

namespace {

// Change this folder name to customize the default CSV destination.
const std::filesystem::path defaultOutputDirectory =
    std::filesystem::path(DEEPNAV_SOURCE_DIR) / "Output clocks";

// =========================== DEMO SETUP ===========================
constexpr double durationDays = 20.0;
constexpr double timeStepSeconds = 1.0;
constexpr double rangeThresholdMeters = 1.0; // Chosen clock-only allocation.
constexpr std::size_t displayStride = 1;
constexpr std::uint64_t localSeed = 42, dsacSeed = 43;
// ================================================================

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Exception = std::invalid_argument, class F> void rejects(F f) {
    try { f(); } catch (const Exception&) { return; }
    throw std::runtime_error("Expected exception was not raised");
}

void historyStorage() {
    static_assert(!std::is_copy_constructible_v<LocalOscillator>);
    static_assert(!std::is_copy_constructible_v<DSAC>);
    static_assert(std::is_move_constructible_v<LocalOscillator>);
    static_assert(std::is_move_constructible_v<DSAC>);
    auto clock = LocalOscillator::fromParameters({});
    require(!clock.historyEnabled() && clock.history() == nullptr, "History must default to disabled");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    clock.recordSample({nan, nan, nan}); // Disabled recording is a no-op.
    clock.clearHistory();
    clock.enableHistory(2, false);
    const auto capacity = clock.history()->samples().capacity();
    const auto* storage = clock.history()->samples().data();
    clock.recordSample({0, 1e-9, 1e-12, nan, nan}); // Unused sigma fields are ignored.
    (void)clock.propagate({}, 1);
    ClockTruthSimulator simulator(clock);
    (void)simulator.step({}, 1);
    require(clock.history()->samples().size() == 1, "Propagation must not record samples");
    rejects([&] { clock.recordSample({0, 0, 0}); });
    rejects([&] { clock.recordSample({-1, 0, 0}); });
    rejects([&] { clock.recordSample({1, nan, 0}); });
    clock.recordSample({1, 2e-9, 1e-12});
    rejects<std::length_error>([&] { clock.recordSample({2, 0, 0}); });
    require(clock.history()->samples().capacity() == capacity
        && clock.history()->samples().data() == storage, "Recording reallocated history");
    clock.clearHistory();
    require(clock.history()->samples().empty() && clock.history()->samples().capacity() == capacity
        && clock.history()->samples().data() == storage, "Clear must preserve allocation");
    clock.recordSample({0, 0, 0});
    auto moved = std::move(clock);
    require(moved.history()->samples().size() == 1 && !clock.historyEnabled(), "Move must transfer history");
    moved.enableHistory(1, true);
    require(moved.history()->samples().empty(), "enableHistory must reset collection");
    rejects([&] { moved.recordSample({0, 0, 0, -1, 0}); });
    rejects([&] { moved.recordSample({0, 0, 0, nan, 0}); });
    moved.disableHistory();
    require(moved.history() == nullptr, "disableHistory must release history");
    moved.enableHistory(0, false);
    rejects<std::length_error>([&] { moved.recordSample({0, 0, 0}); });
    require(ClockHistory::requiredCapacity(10, 3) == 5, "Stride capacity incorrect");
    require(ClockHistory::requiredCapacity(9, 3) == 4, "Exact stride capacity incorrect");
    require(ClockHistory::requiredCapacity(0, 3) == 1, "Initial-only capacity incorrect");
    rejects([&] { (void)ClockHistory::requiredCapacity(1, 0); });
    rejects<std::overflow_error>([&] {
        (void)ClockHistory::requiredCapacity(std::numeric_limits<std::size_t>::max(), 1);
    });
    std::cout << "PASS: optional/movable history, validation, limits and allocation preservation\n";
}

void recordingRuns() {
    auto clock = DSAC::fromParameters({0, 1e-24, 3e-26});
    ClockTruthSimulator simulator(clock, 42);
    clock.enableHistory(5, true);
    rejects([&] { (void)simulator.run(10, 2, {}, 3); });
    require(clock.history()->samples().empty(), "Failed preparation must not record");
    const ClockCovariance p0 = ClockCovariance::Zero();
    const auto final = simulator.run(10, 2, {}, 3, p0);
    const auto& samples = clock.history()->samples();
    require(samples.size() == 5, "Stride must include initial and final samples exactly once");
    const double times[]{0, 6, 12, 18, 20};
    for (std::size_t i = 0; i < samples.size(); ++i) {
        require(samples[i].time_s == times[i], "Incorrect recorded timestamp");
        const double t = times[i];
        const double expected = std::sqrt(1e-24 * t + 3e-26 * t * t * t / 3);
        require(std::abs(samples[i].sigma_bias_s - expected) <= expected * 1e-12, "Incorrect bias sigma");
        require(std::abs(samples[i].sigma_fractional_frequency - std::sqrt(3e-26 * t)) < 1e-25,
            "Incorrect frequency sigma");
    }
    require(samples.back().bias_s == final.bias_s, "Final truth sample mismatch");
    const auto reference = DSAC::fromParameters(clock.parameters());
    ClockTruthSimulator unrecorded(reference, 42);
    ClockState manual;
    for (int i = 0; i < 10; ++i) manual = unrecorded.step(manual, 2);
    require(manual.bias_s == final.bias_s && manual.fractional_frequency == final.fractional_frequency,
        "Recording changed the truth realization");
    rejects([&] { (void)simulator.run(10, 2, {}, 3, p0); });
    clock.enableHistory(4, false);
    rejects<std::length_error>([&] { (void)simulator.run(10, 2, {}, 3); });
    require(clock.history()->samples().empty(), "Insufficient capacity must fail before recording");
    (void)simulator.run(9, 2, {}, 3);
    require(clock.history()->samples().size() == 4, "Final multiple of stride was duplicated");
    clock.enableHistory(3, false);
    const double variable[]{0.25, 1.5, 0.125};
    (void)simulator.run(variable, {}, 2);
    require(clock.history()->samples()[1].time_s == 1.75
        && clock.history()->samples().back().time_s == 1.875, "Variable dt timestamps were reconstructed");
    clock.enableHistory(1, false);
    (void)simulator.run(0, 1);
    require(clock.history()->samples().size() == 1, "Zero steps must retain initial sample");
    clock.clearHistory();
    const Clocks& const_clock = clock;
    ClockTruthSimulator read_only(const_clock);
    rejects([&] { (void)read_only.run(0, 1); });
    clock.disableHistory();
    (void)read_only.run(2, 1);
    rejects([&] { (void)read_only.run(1, 1, {}, 0); });
    std::cout << "PASS: accepted-step recording, stride endpoints, P0 requirements and variable elapsed time\n";
}

void csvExport(const std::filesystem::path& directory) {
    auto clock = LocalOscillator::fromParameters({});
    clock.enableHistory(2, true);
    const double bias = std::nextafter(1e-9, 2e-9);
    clock.recordSample({0, bias, -1e-12, 2e-10, 3e-15});
    clock.recordSample({2, 2e-9, -2e-12, 3e-10, 4e-15});
    const auto path = directory / "precision.csv";
    writeClockCsv(path, *clock.history(), ClockCsvWriteMode::Overwrite);
    rejects<std::runtime_error>([&] { writeClockCsv(path, *clock.history()); });
    std::ifstream input(path);
    std::string header, row;
    std::getline(input, header);
    require(header == "time_s,bias_s,fractional_frequency,sigma_bias_s,sigma_fractional_frequency",
        "Wrong covariance CSV header");
    std::getline(input, row);
    std::istringstream fields(row);
    std::string field;
    std::getline(fields, field, ',');
    require(std::stod(field) == 0, "Time units changed on export");
    std::getline(fields, field, ',');
    require(std::stod(field) == bias, "CSV did not round-trip double precision");
    std::getline(fields, field, ',');
    require(std::stod(field) == -1e-12, "Frequency units changed on export");
    clock.enableHistory(0, false);
    writeClockCsv(path, *clock.history(), ClockCsvWriteMode::Overwrite);
    std::ifstream empty(path);
    std::getline(empty, header);
    require(header == "time_s,bias_s,fractional_frequency", "Wrong no-covariance header");
    require(!std::getline(empty, row), "Overwrite appended old samples");
    rejects<std::runtime_error>([&] { writeClockCsv(directory / "missing" / "file.csv", *clock.history()); });
    rejects<std::runtime_error>([&] { writeClockCsv(directory, *clock.history(), ClockCsvWriteMode::Overwrite); });
    if (std::filesystem::exists("/dev/full"))
        rejects<std::runtime_error>([&] { writeClockCsv("/dev/full", *clock.history(), ClockCsvWriteMode::Overwrite); });
    std::cout << "PASS: exact CSV schema, precision, empty export, overwrite policy and I/O failures\n";
}

void exportDeterministicMean(const Clocks& clock, ClockState state,
    const std::filesystem::path& path, std::size_t steps, double dt, std::size_t stride) {
    auto mean = LocalOscillator::fromParameters(clock.parameters());
    mean.enableHistory(ClockHistory::requiredCapacity(steps, stride), false);
    mean.recordSample({0, state.bias_s, state.fractional_frequency});
    double elapsed = 0;
    for (std::size_t i = 1; i <= steps; ++i) {
        state = mean.propagate(state, dt); // Expected state; no random increments.
        elapsed += dt;
        if (i % stride == 0 || i == steps)
            mean.recordSample({elapsed, state.bias_s, state.fractional_frequency});
    }
    writeClockCsv(path, *mean.history(), ClockCsvWriteMode::Overwrite);
}

void exportPlotFixtures(const std::filesystem::path& directory) {
    auto local = LocalOscillator::representativeUsoWhiteFmWithAging();
    auto dsac = DSAC::dayMatchedWhiteFmBaseline();
    const double durationSeconds = durationDays * 86400.0;
    const double requestedSteps = durationSeconds / timeStepSeconds;
    if (!std::isfinite(durationSeconds) || durationSeconds <= 0 || !std::isfinite(timeStepSeconds)
        || timeStepSeconds <= 0 || !std::isfinite(requestedSteps)
        || requestedSteps >= static_cast<double>(std::numeric_limits<std::size_t>::max())
        || std::abs(requestedSteps - std::round(requestedSteps)) > 1e-9
        || displayStride == 0 || !std::isfinite(rangeThresholdMeters) || rangeThresholdMeters <= 0)
        throw std::invalid_argument("Demo setup requires positive duration/dt/threshold, integral steps and stride >= 1.");
    const std::size_t steps = static_cast<std::size_t>(std::round(requestedSteps)), stride = displayStride;
    const double dt = timeStepSeconds;
    // Ideal initial bias/frequency calibration. Customize residual state and
    // P0 here independently of the oscillator's diffusion intensities.
    const ClockState localInitial{}, dsacInitial{};
    const ClockCovariance localP0 = ClockCovariance::Zero();
    const ClockCovariance dsacP0 = ClockCovariance::Zero();
    // Keep every step for Allan; display CSVs are selected from the same truth.
    local.enableHistory(ClockHistory::requiredCapacity(steps, 1), true);
    dsac.enableHistory(ClockHistory::requiredCapacity(steps, 1), true);
    ClockTruthSimulator a(local, localSeed), b(dsac, dsacSeed);
    (void)a.run(steps, dt, localInitial, 1, localP0);
    (void)b.run(steps, dt, dsacInitial, 1, dsacP0);
    writeClockCsv(directory / "local_allan.csv", *local.history(), ClockCsvWriteMode::Overwrite);
    writeClockCsv(directory / "dsac_allan.csv", *dsac.history(), ClockCsvWriteMode::Overwrite);
    for (const auto& entry : {std::pair<const Clocks*, const char*>{&local, "local.csv"},
                              std::pair<const Clocks*, const char*>{&dsac, "dsac.csv"}}) {
        auto display = LocalOscillator::fromParameters(entry.first->parameters());
        display.enableHistory(ClockHistory::requiredCapacity(steps, stride), true);
        const auto& samples = entry.first->history()->samples();
        for (std::size_t i = 0; i < samples.size(); ++i)
            if (i % stride == 0 || i + 1 == samples.size()) display.recordSample(samples[i]);
        writeClockCsv(directory / entry.second, *display.history(), ClockCsvWriteMode::Overwrite);
    }
    exportDeterministicMean(local, localInitial, directory / "local_mean.csv", steps, dt, stride);
    exportDeterministicMean(dsac, dsacInitial, directory / "dsac_mean.csv", steps, dt, stride);
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "Isolated clock-only comparison: ideal ground clock; no measurement noise.\n"
              << "  " << local.configurationName() << ": D=" << local.parameters().frequency_drift_per_s
              << " /s, q_b=" << local.parameters().q_bias_s << " s, q_y=0 /s, seed=" << localSeed
              << "; aging plus simplified white FM, not a full device spectrum fit.\n"
              << "  " << dsac.configurationName() << ": D=" << dsac.parameters().frequency_drift_per_s
              << " /s, q_b=" << dsac.parameters().q_bias_s
              << " s, q_y=0 /s, seed=" << dsacSeed << "; effective approximation, not a fitted DSAC hardware model.\n"
              << "  local initial: b0=" << localInitial.bias_s << " s, y0=" << localInitial.fractional_frequency
              << "; P0=[" << localP0(0,0) << ',' << localP0(0,1) << ';' << localP0(1,0) << ',' << localP0(1,1) << "]\n"
              << "  DSAC initial: b0=" << dsacInitial.bias_s << " s, y0=" << dsacInitial.fractional_frequency
              << "; P0=[" << dsacP0(0,0) << ',' << dsacP0(0,1) << ';' << dsacP0(1,0) << ',' << dsacP0(1,1) << "]\n"
              << "  Sampling: dt=" << dt << " s, duration=" << steps * dt
              << " s, display stride=" << stride << ", Allan stride=1\n"
              << "  Allocation: " << rangeThresholdMeters << " m of clock-only range error, not a navigation requirement.\n"
              << "  Ideal initial calibration; ideal ground reference; constant aging uncompensated; simplified white FM; flicker/other long-term noise omitted.\n"
              << "  Simplified holdover comparison; sampled crossings are not mandatory ground-contact intervals.\n"
              << "  Sources: https://doi.org/10.1029/2025RS008244; https://doi.org/10.1038/s41586-021-03571-7\n"
              << "  One-day ADEV reference: Tjoelker (2021), DSAC results, slide 21.\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc > 2) {
        std::cerr << "Usage: test_clock_history [output-directory]\n";
        return 2;
    }
    const auto directory = argc == 2 ? std::filesystem::path(argv[1]) : defaultOutputDirectory;
    const auto scratch = std::filesystem::temp_directory_path() / ("deepnav-clocks-"
        + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        std::filesystem::create_directories(scratch);
        historyStorage();
        recordingRuns();
        csvExport(scratch);
        std::filesystem::remove_all(scratch);
        std::filesystem::create_directories(directory);
        exportPlotFixtures(directory);
        std::cout << "Saved reference-model clock CSVs to " << std::filesystem::absolute(directory) << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::error_code cleanupError;
        std::filesystem::remove_all(scratch, cleanupError);
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
