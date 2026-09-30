#pragma once

#include "dynamics/CartesianState.hpp"
#include "dynamics/StateProvider.hpp"
#include "opnav/Types.hpp"

#include <cstddef>

namespace fd::opnav {

struct LightTimeSolverOptions {
    double toleranceSeconds{1.0e-9};
    std::size_t maxIterations{8};
};

struct LightTimeSolution {
    TdbEpoch receptionEpoch;
    TdbEpoch emissionEpoch;

    double lightTimeSeconds;
    double residualSeconds;

    dynamics::CartesianState targetStateAtEmission;

    // T(t) = S(t - tau) - R(t)
    Eigen::Vector3d geometricLineOfSightKm;

    std::size_t iterations;
};

class LightTimeSolver {
public:
    explicit constexpr LightTimeSolver(
        LightTimeSolverOptions options = {})
        : options_(options) {}

    [[nodiscard]]
    LightTimeSolution solve(
        TdbEpoch receptionEpoch,
        const dynamics::CartesianState& observerStateAtReception,
        const dynamics::StateProvider& targetProvider) const;

private:
    static constexpr double speedOfLightKmPerSecond =
        299792.458;

    LightTimeSolverOptions options_;
};

} // namespace fd::opnav