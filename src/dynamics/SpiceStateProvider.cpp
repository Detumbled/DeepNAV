#include "dynamics/SpiceStateProvider.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <cmath>
#include <utility>

namespace fd::dynamics {
SpiceStateProvider::SpiceStateProvider(std::string target) : target_(std::move(target)) {
    if (target_.empty()) throw std::invalid_argument("SPICE target cannot be empty.");
}
CartesianState SpiceStateProvider::stateAt(opnav::TdbEpoch epoch) const {
    if (!std::isfinite(epoch.secondsPastJ2000))
        throw std::invalid_argument("SPICE epoch must be finite.");
    od::SpiceErrorModeGuard guard;
    SpiceDouble state[6], lightTime;
    spkezr_c(target_.c_str(), epoch.secondsPastJ2000, "J2000", "NONE", "SSB", state, &lightTime);
    od::throwIfSpiceFailed("Cannot read geometric state of " + target_);
    return {Eigen::Map<Eigen::Vector3d>(state), Eigen::Map<Eigen::Vector3d>(state + 3)};
}
} // namespace fd::dynamics
