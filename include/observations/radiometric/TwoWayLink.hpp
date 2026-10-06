#pragma once

#include "dynamics/StateProvider.hpp"
#include "opnav/core/LightTimeSolver.hpp"
#include <functional>

namespace fd::observations::radiometric {
inline constexpr double speedOfLightKmPerSecond = 299792.458;

struct LegSolution {
    double receptionTdb{}, emissionTdb{};
    double geometricSeconds{}, correctionSeconds{}, residualSeconds{};
    dynamics::CartesianState emitter, receiver;
    std::size_t iterations{};
    [[nodiscard]] double durationSeconds() const { return geometricSeconds + correctionSeconds; }
};
using LegSolver = std::function<LegSolution(double, const dynamics::StateProvider&,
                                           const dynamics::StateProvider&)>;
// Delay in coordinate seconds at the endpoints. Supply group delay for range,
// phase delay for Doppler; dispersive media generally require different models.
using DelayModel = std::function<double(const LegSolution&, const dynamics::StateProvider&,
                                        const dynamics::StateProvider&)>;
[[nodiscard]] LegSolver libraryReceptionSolver(opnav::LightTimeSolverOptions options = {},
                                               DelayModel delay = {});

struct LinkConfig {
    // Must be supplied by the mission/transponder configuration.
    double turnaroundRatio{};
    // Effective constant non-dispersive delay in TDB seconds, NOT a flight calibration.
    double transponderDelaySeconds{0};
};
struct TwoWaySolution {
    LegSolution uplink, downlink;
    double spacecraftReceptionTdb{}, spacecraftTransmissionTdb{};
    double transponderDelaySeconds{};
    [[nodiscard]] double roundTripSeconds() const;
    // Conventional c * RTLT / 2; this is not either leg's instantaneous distance.
    [[nodiscard]] double rangeKm() const;
};
struct FrequencyRamp {
    double referenceTdb{}, frequencyHz{}, rateHzPerSecond{};
};
// Returns TDB minus station clock reading (e.g. TT), in seconds. An empty callback
// means coordinate TDB frequencies. Local clock errors can be added here.
using ClockOffset = std::function<double(double)>;
struct DopplerCount {
    TwoWaySolution start, end;
    double receiveIntervalSeconds{}, transmitIntervalSeconds{};
    double receivedCycles{}, meanReceivedHz{}, residualHz{};
};

class TwoWayLink {
  public:
    TwoWayLink(LegSolver solver, LinkConfig config);
    [[nodiscard]] TwoWaySolution solve(double receptionTdb,
                                       const dynamics::StateProvider& station,
                                       const dynamics::StateProvider& spacecraft) const;
    // End-tagged count; positive residual means received frequency above reference.
    // Spacecraft oscillator noise cancels in this ideal coherently locked model.
    [[nodiscard]] DopplerCount count(double receptionTdb, double intervalTdbSeconds,
                                     const dynamics::StateProvider& station,
                                     const dynamics::StateProvider& spacecraft,
                                     const FrequencyRamp& uplink, double referenceHz,
                                     ClockOffset groundClock = {}) const;
  private:
    LegSolver solver_;
    LinkConfig config_;
};
} // namespace fd::observations::radiometric
