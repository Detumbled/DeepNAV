#include "observations/radiometric/OneWayLink.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::observations::radiometric {
OneWayLink::OneWayLink(LegSolver group, LegSolver phase)
    : group_(std::move(group)), phase_(phase ? std::move(phase) : group_) {
    if (!group_)
        throw std::invalid_argument("One-way leg solver required.");
}
LegSolution OneWayLink::solve(double t, const dynamics::StateProvider &tx,
                              const dynamics::StateProvider &rx) const {
    if (!std::isfinite(t))
        throw std::invalid_argument("Reception epoch must be finite.");
    return group_(t, tx, rx);
}
OneWayCount OneWayLink::count(long double tag, double interval, const dynamics::StateProvider &tx,
                              const dynamics::StateProvider &rx, const EndpointClock &ground,
                              const EndpointClock &onboard, const FrequencyRamp &ramp,
                              double referenceHz) const {
    if (!std::isfinite(interval) || interval <= 0 || !std::isfinite(ramp.referenceTdb) ||
        !std::isfinite(ramp.frequencyHz) || ramp.frequencyHz <= 0 ||
        !std::isfinite(ramp.rateHzPerSecond) || !std::isfinite(referenceHz) || referenceHz <= 0)
        throw std::invalid_argument("Invalid one-way count or frequency ramp.");
    const double end = onboard.coordinateEpoch(tag);
    // Solve a SHORT interval rather than subtracting absolute clock epochs.
    // Keep the unquantized interval for GHz phase integration; trajectory query
    // epochs alone are rounded to the StateProvider's double TDB representation.
    long double dt = interval;
    double begin{};
    bool converged = false;
    for (int i = 0; i < 32; ++i) {
        begin = static_cast<double>(static_cast<long double>(end) - dt);
        const long double next = interval - (static_cast<long double>(onboard.correction(end)) -
                                             onboard.correction(begin));
        if (!std::isfinite(next) || next <= 0)
            throw std::runtime_error("Non-positive receiver coordinate interval.");
        if (std::abs(next - dt) <= 1e-12L) {
            dt = next;
            converged = true;
            break;
        }
        dt = next;
    }
    if (!converged)
        throw std::runtime_error("Receiver count interval did not converge.");
    begin = static_cast<double>(static_cast<long double>(end) - dt);
    if (begin == end)
        throw std::invalid_argument("Count interval below epoch resolution.");
    OneWayCount result;
    result.start = phase_(begin, tx, rx);
    result.end = phase_(end, tx, rx);
    result.rangeEnd = group_(end, tx, rx);
    // Compute the small interval difference directly, including receiver-clock
    // error, rather than subtracting two ~60 s intervals (long double is only
    // double precision on some platforms).
    const long double intervalDifference =
        -(static_cast<long double>(onboard.correction(end)) - onboard.correction(begin)) -
        (static_cast<long double>(result.end.durationSeconds()) - result.start.durationSeconds()) +
        (static_cast<long double>(ground.correction(result.end.emissionTdb)) -
         ground.correction(result.start.emissionTdb));
    const long double txInterval = interval + intervalDifference;
    if (!std::isfinite(txInterval) || txInterval <= 0)
        throw std::runtime_error("Non-positive transmitter clock interval.");
    const long double txStart =
        (static_cast<long double>(end) - ramp.referenceTdb) - dt - result.start.durationSeconds() +
        ground.correction(result.start.emissionTdb) - ground.correction(ramp.referenceTdb);
    const long double meanCarrier =
        ramp.frequencyHz + ramp.rateHzPerSecond * (txStart + .5L * txInterval);
    if (!std::isfinite(meanCarrier) ||
        meanCarrier - .5L * std::abs(ramp.rateHzPerSecond) * txInterval <= 0)
        throw std::runtime_error("Ramp leaves positive frequency domain.");
    result.receiveCoordinateIntervalSeconds = dt;
    result.receiveClockIntervalSeconds = interval;
    result.transmitClockIntervalSeconds = txInterval;
    result.receivedCycles = meanCarrier * txInterval;
    result.residualHz = (meanCarrier - referenceHz) + meanCarrier * (intervalDifference / interval);
    result.pseudorangeKm =
        speedOfLightKmPerSecond * (result.rangeEnd.durationSeconds() + onboard.correction(end) -
                                   ground.correction(result.rangeEnd.emissionTdb));
    if (!std::isfinite(result.residualHz) || !std::isfinite(result.pseudorangeKm) ||
        !std::isfinite(result.receivedCycles))
        throw std::runtime_error("One-way observation overflowed.");
    return result;
}
} // namespace fd::observations::radiometric
