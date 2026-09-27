#include "dynamics/LinearStateProvider.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::dynamics {

LinearStateProvider::LinearStateProvider(
    opnav::TdbEpoch referenceEpoch,
    CartesianState referenceState,
    std::string referenceFrame,
    std::string referenceOrigin)
    : referenceEpoch_(referenceEpoch),
      referenceState_(std::move(referenceState)),
      referenceFrame_(std::move(referenceFrame)),
      referenceOrigin_(std::move(referenceOrigin)) {

    if (!std::isfinite(referenceEpoch_.secondsPastJ2000)) {
        throw std::invalid_argument(
            "Reference epoch must be finite");
    }

    if (!referenceState_.allFinite()) {
        throw std::invalid_argument(
            "Reference state must contain finite values");
    }

    if (referenceFrame_.empty() || referenceOrigin_.empty()) {
        throw std::invalid_argument(
            "Reference frame and origin cannot be empty");
    }
}

CartesianState LinearStateProvider::stateAt(
    opnav::TdbEpoch epoch) const {

    if (!std::isfinite(epoch.secondsPastJ2000)) {
        throw std::invalid_argument(
            "Requested epoch must be finite");
    }

    const double deltaTimeSeconds =
        epoch.secondsPastJ2000 -
        referenceEpoch_.secondsPastJ2000;

    CartesianState result;
    result.positionKm =
        referenceState_.positionKm +
        referenceState_.velocityKmPerSec * deltaTimeSeconds;
    result.velocityKmPerSec =
        referenceState_.velocityKmPerSec;

    if (!result.allFinite()) {
        throw std::runtime_error(
            "Linear propagation produced a non-finite state");
    }

    return result;
}

std::string_view
LinearStateProvider::referenceFrame() const noexcept {
    return referenceFrame_;
}

std::string_view
LinearStateProvider::referenceOrigin() const noexcept {
    return referenceOrigin_;
}

} // namespace fd::dynamics
