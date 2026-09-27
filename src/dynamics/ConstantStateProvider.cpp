#include "dynamics/ConstantStateProvider.hpp"

#include <stdexcept>
#include <utility>

namespace fd::dynamics {

ConstantStateProvider::ConstantStateProvider(
    CartesianState state,
    std::string referenceFrame,
    std::string referenceOrigin)
    : state_(std::move(state)),
      referenceFrame_(std::move(referenceFrame)),
      referenceOrigin_(std::move(referenceOrigin)) {

    if (!state_.allFinite()) {
        throw std::invalid_argument(
            "Constant state must contain finite values");
    }

    if (referenceFrame_.empty() || referenceOrigin_.empty()) {
        throw std::invalid_argument(
            "Reference frame and origin cannot be empty");
    }
}

CartesianState ConstantStateProvider::stateAt(
    opnav::TdbEpoch epoch) const {
    (void)epoch;
    return state_;
}

std::string_view
ConstantStateProvider::referenceFrame() const noexcept {
    return referenceFrame_;
}

std::string_view
ConstantStateProvider::referenceOrigin() const noexcept {
    return referenceOrigin_;
}

} // namespace fd::dynamics
