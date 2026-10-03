#pragma once

#include "opnav/core/LightTimeSolver.hpp"

namespace fd::opnav::core {

struct ApparentDirection {
    LightTimeSolution lightTime;
    Eigen::Vector3d vectorKm;
    Eigen::Vector3d unitDirection;
};

// Both states must be geometric (no corrections already applied), in the same
// barycentric inertial frame. Observer velocity must be relative to the SSB.
// Owen, Section 4.2: A(t) = S(t - tau) - R(t) + tau * v_observer(t).
[[nodiscard]] ApparentDirection computeApparentDirection(
    TdbEpoch receptionEpoch,
    const dynamics::CartesianState& observerStateAtReception,
    const dynamics::StateProvider& targetProvider,
    LightTimeSolverOptions options = {});

} // namespace fd::opnav::core
