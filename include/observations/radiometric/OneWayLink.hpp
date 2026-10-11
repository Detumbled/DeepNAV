#pragma once
#include "observations/radiometric/EndpointClock.hpp"

namespace fd::observations::radiometric {
struct OneWayCount {
    LegSolution start, end, rangeEnd;
    double receiveCoordinateIntervalSeconds{}, receiveClockIntervalSeconds{};
    double transmitClockIntervalSeconds{}, receivedCycles{}, residualHz{}, pseudorangeKm{};
};
class OneWayLink {
  public:
    // Separate group and phase solvers support dispersive corrections. If phase
    // is empty, group is reused (non-dispersive link).
    explicit OneWayLink(LegSolver group, LegSolver phase = {});
    [[nodiscard]] LegSolution solve(double receptionTdb, const dynamics::StateProvider &tx,
                                    const dynamics::StateProvider &rx) const;
    // End tag and count duration are RECEIVER clock seconds. The ramp is Hz per
    // transmitter-clock second, anchored at ramp.referenceTdb. Positive residual
    // means received phase rate exceeds the receiver's local referenceHz.
    [[nodiscard]] OneWayCount count(long double receptionClockSeconds, double countClockSeconds,
                                    const dynamics::StateProvider &tx,
                                    const dynamics::StateProvider &rx,
                                    const EndpointClock &transmitterClock,
                                    const EndpointClock &receiverClock, const FrequencyRamp &ramp,
                                    double referenceHz) const;

  private:
    LegSolver group_, phase_;
};
} // namespace fd::observations::radiometric
