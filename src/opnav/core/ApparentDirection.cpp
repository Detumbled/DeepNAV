#include "opnav/core/ApparentDirection.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::opnav::core {

ApparentDirection computeApparentDirection(
    TdbEpoch receptionEpoch,
    const dynamics::CartesianState& observerStateAtReception,
    const dynamics::StateProvider& targetProvider,
    LightTimeSolverOptions options) {
    auto solution = LightTimeSolver(options).solve(
        receptionEpoch, observerStateAtReception, targetProvider);
    const Eigen::Vector3d apparent = solution.geometricLineOfSightKm
        + solution.lightTimeSeconds * observerStateAtReception.velocityKmPerSec;
    const double length = apparent.stableNorm();
    if (!(length > 0.0) || !std::isfinite(length)) {
        throw std::runtime_error("Apparent direction must be finite and nonzero.");
    }
    return {std::move(solution), apparent, apparent / length};
}

} // namespace fd::opnav::core
