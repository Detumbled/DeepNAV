#include "observations/radiometric/SpiceLightTime.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <algorithm>
#include <cmath>

namespace fd::observations::radiometric {
LegSolution spiceReceptionLeg(double reception, const dynamics::StateProvider &emitter,
                              const dynamics::StateProvider &receiver) {
    const auto *tx = dynamic_cast<const dynamics::SpiceStateProvider *>(&emitter);
    const auto *rx = dynamic_cast<const dynamics::SpiceStateProvider *>(&receiver);
    if (!tx || !rx || !std::isfinite(reception))
        throw std::invalid_argument("CSPICE CN requires finite epoch and SPICE endpoints.");
    od::SpiceErrorModeGuard guard;
    SpiceDouble state[6], lt;
    spkezr_c(tx->target().c_str(), reception, "J2000", "CN", rx->target().c_str(), state, &lt);
    od::throwIfSpiceFailed("CSPICE CN light time");
    LegSolution leg;
    leg.receptionTdb = reception;
    leg.emissionTdb = reception - lt;
    leg.emitter = emitter.stateAt(opnav::TdbEpoch{leg.emissionTdb});
    leg.receiver = receiver.stateAt(opnav::TdbEpoch{reception});
    leg.geometricSeconds = lt;
    leg.residualSeconds = std::abs(lt - (leg.receiver.positionKm - leg.emitter.positionKm).norm() /
                                            speedOfLightKmPerSecond);
    return leg;
}
double spiceStationElevationRadians(const std::string &station, double epoch,
                                    const Eigen::Vector3d &direction) {
    if (station.empty() || !std::isfinite(epoch) || !direction.allFinite() || direction.norm() == 0)
        throw std::invalid_argument("Invalid station elevation geometry.");
    od::SpiceErrorModeGuard guard;
    SpiceDouble rotation[3][3];
    pxform_c("J2000", (station + "_TOPO").c_str(), epoch, rotation);
    od::throwIfSpiceFailed("Station topocentric transform");
    const Eigen::Vector3d local =
        Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(&rotation[0][0]) * direction;
    return std::asin(std::clamp(local.z() / local.norm(), -1.0, 1.0));
}
SpiceTdbClock::SpiceTdbClock() {
    od::SpiceErrorModeGuard guard;
    const auto constant = [](const char *name, int count, double *values) {
        SpiceInt n;
        SpiceBoolean found;
        gdpool_c(name, 0, count, &n, values, &found);
        od::throwIfSpiceFailed("Read LSK DELTET constants");
        if (!found || n != count)
            throw std::runtime_error("Incomplete LSK DELTET model.");
    };
    double m[2];
    constant("DELTET/K", 1, &k_);
    constant("DELTET/EB", 1, &eb_);
    constant("DELTET/M", 2, m);
    meanAnomaly_ = m[0];
    anomalyRate_ = m[1];
    if (!std::isfinite(k_) || !std::isfinite(eb_) || !std::isfinite(meanAnomaly_) ||
        !std::isfinite(anomalyRate_))
        throw std::runtime_error("Non-finite LSK DELTET constants.");
}
double SpiceTdbClock::operator()(double epoch) const {
    if (!std::isfinite(epoch))
        throw std::invalid_argument("Clock epoch must be finite.");
    const double anomaly = meanAnomaly_ + anomalyRate_ * epoch;
    return k_ * std::sin(anomaly + eb_ * std::sin(anomaly));
}
double spiceTdbMinusTt(double epoch) { return SpiceTdbClock{}(epoch); }
} // namespace fd::observations::radiometric
