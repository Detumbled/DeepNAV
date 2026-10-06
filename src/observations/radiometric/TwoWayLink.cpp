#include "observations/radiometric/TwoWayLink.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::observations::radiometric {
namespace {
void commonFrame(const dynamics::StateProvider &a, const dynamics::StateProvider &b) {
    if (a.referenceFrame() != b.referenceFrame() || a.referenceOrigin() != b.referenceOrigin())
        throw std::invalid_argument(
            "Light-time endpoints must share an inertial frame and origin.");
}
LegSolution endpoints(double reception, double duration, const dynamics::StateProvider &tx,
                      const dynamics::CartesianState &receiver) {
    LegSolution leg;
    leg.receptionTdb = reception;
    leg.emissionTdb = reception - duration;
    leg.emitter = tx.stateAt(opnav::TdbEpoch{leg.emissionTdb});
    leg.receiver = receiver;
    if (!leg.emitter.allFinite() || !leg.receiver.allFinite())
        throw std::runtime_error("Non-finite radio endpoint state.");
    leg.geometricSeconds =
        (leg.receiver.positionKm - leg.emitter.positionKm).norm() / speedOfLightKmPerSecond;
    return leg;
}
} // namespace

LegSolver libraryReceptionSolver(opnav::LightTimeSolverOptions options, DelayModel delay) {
    return [options, delay = std::move(delay)](double reception, const dynamics::StateProvider &tx,
                                               const dynamics::StateProvider &rx) {
        commonFrame(tx, rx);
        // Reception is fixed throughout this leg. Keep the receiver state and
        // reuse the geometric solver's converged emitter instead of querying again.
        const auto receiver = rx.stateAt(opnav::TdbEpoch{reception});
        const auto seed =
            opnav::LightTimeSolver(options).solve(opnav::TdbEpoch{reception}, receiver, tx);
        LegSolution leg;
        leg.receptionTdb = reception;
        leg.emissionTdb = seed.emissionEpoch.secondsPastJ2000;
        leg.emitter = seed.targetStateAtEmission;
        leg.receiver = receiver;
        leg.geometricSeconds = seed.geometricLineOfSightKm.norm() / speedOfLightKmPerSecond;
        leg.residualSeconds = seed.residualSeconds;
        leg.iterations = seed.iterations;
        if (!delay)
            return leg;
        double duration = seed.lightTimeSeconds;
        for (std::size_t i = 1; i <= options.maxIterations; ++i) {
            leg.correctionSeconds = delay(leg, tx, rx);
            const double next = leg.durationSeconds();
            if (!std::isfinite(next) || next <= 0)
                throw std::runtime_error("Invalid corrected light time.");
            leg.residualSeconds = std::abs(next - duration);
            leg.iterations = seed.iterations + i;
            if (leg.residualSeconds <= options.toleranceSeconds)
                return leg;
            duration = next;
            if (i < options.maxIterations)
                leg = endpoints(reception, duration, tx, receiver);
        }
        throw std::runtime_error("Corrected light-time iteration did not converge.");
    };
}

double TwoWaySolution::roundTripSeconds() const {
    return uplink.durationSeconds() + transponderDelaySeconds + downlink.durationSeconds();
}
double TwoWaySolution::rangeKm() const { return .5 * speedOfLightKmPerSecond * roundTripSeconds(); }
TwoWayLink::TwoWayLink(LegSolver solver, LinkConfig config)
    : solver_(std::move(solver)), config_(config) {
    if (!solver_ || !std::isfinite(config.turnaroundRatio) || config.turnaroundRatio <= 0 ||
        !std::isfinite(config.transponderDelaySeconds) || config.transponderDelaySeconds < 0)
        throw std::invalid_argument("Invalid coherent two-way link configuration.");
}
TwoWaySolution TwoWayLink::solve(double reception, const dynamics::StateProvider &station,
                                 const dynamics::StateProvider &spacecraft) const {
    if (!std::isfinite(reception))
        throw std::invalid_argument("Reception epoch must be finite.");
    TwoWaySolution result;
    result.downlink = solver_(reception, spacecraft, station);
    result.spacecraftTransmissionTdb = result.downlink.emissionTdb;
    result.spacecraftReceptionTdb =
        result.spacecraftTransmissionTdb - config_.transponderDelaySeconds;
    // On the uplink Cassini is the receiver, and the station is the retarded emitter.
    result.uplink = solver_(result.spacecraftReceptionTdb, station, spacecraft);
    result.transponderDelaySeconds = config_.transponderDelaySeconds;
    return result;
}
DopplerCount TwoWayLink::count(double reception, double interval,
                               const dynamics::StateProvider &station,
                               const dynamics::StateProvider &spacecraft, const FrequencyRamp &ramp,
                               double referenceHz, ClockOffset clock) const {
    if (!std::isfinite(interval) || interval <= 0 || !std::isfinite(ramp.referenceTdb) ||
        !std::isfinite(ramp.frequencyHz) || ramp.frequencyHz <= 0 ||
        !std::isfinite(ramp.rateHzPerSecond) || !std::isfinite(referenceHz) || referenceHz <= 0)
        throw std::invalid_argument("Invalid Doppler count or frequency ramp.");
    const auto offset = [&](double t) {
        const double value = clock ? clock(t) : 0;
        if (!std::isfinite(value))
            throw std::runtime_error("Non-finite station clock offset.");
        return value;
    };
    if (!std::isfinite(reception) || !std::isfinite(reception - interval) ||
        reception - interval == reception)
        throw std::invalid_argument("Count interval is not representable at this TDB epoch.");
    // Use the actual representable endpoints, including for fractional intervals.
    interval = reception - (reception - interval);
    DopplerCount result;
    result.start = solve(reception - interval, station, spacecraft);
    result.end = solve(reception, station, spacecraft);
    // Difference the short delays, not large absolute emission epochs: otherwise
    // epoch quantization is amplified by the GHz carrier into spurious Doppler.
    const double txStartOffset = offset(result.start.uplink.emissionTdb);
    const double txEndOffset = offset(result.end.uplink.emissionTdb);
    const long double txInterval =
        static_cast<long double>(interval) -
        (result.end.roundTripSeconds() - result.start.roundTripSeconds()) -
        (txEndOffset - txStartOffset);
    const long double rxInterval =
        static_cast<long double>(interval) - (offset(reception) - offset(reception - interval));
    if (txInterval <= 0 || rxInterval <= 0)
        throw std::runtime_error("Non-positive station clock count interval.");
    const long double txStart = static_cast<long double>(reception - interval) - ramp.referenceTdb -
                                result.start.roundTripSeconds() -
                                (txStartOffset - offset(ramp.referenceTdb));
    const long double meanUplink =
        ramp.frequencyHz + ramp.rateHzPerSecond * (txStart + .5L * txInterval);
    if (meanUplink <= 0 || !std::isfinite(meanUplink))
        throw std::runtime_error("Frequency ramp leaves the positive carrier domain.");
    const long double carrier = config_.turnaroundRatio * meanUplink;
    const long double residual =
        (carrier - referenceHz) + carrier * ((txInterval - rxInterval) / rxInterval);
    result.receiveIntervalSeconds = rxInterval;
    result.transmitIntervalSeconds = txInterval;
    result.receivedCycles = carrier * txInterval;
    result.residualHz = residual;
    result.meanReceivedHz = referenceHz + residual;
    return result;
}
} // namespace fd::observations::radiometric
