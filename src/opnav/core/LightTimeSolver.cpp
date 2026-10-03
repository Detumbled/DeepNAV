#include "opnav/core/LightTimeSolver.hpp"

#include <cmath>
#include <stdexcept>

namespace fd::opnav {

LightTimeSolution LightTimeSolver::solve(
    TdbEpoch receptionEpoch,
    const dynamics::CartesianState& observerStateAtReception,
    const dynamics::StateProvider& targetProvider) const {

    if (!std::isfinite(receptionEpoch.secondsPastJ2000)
        || !observerStateAtReception.allFinite()) {
        throw std::invalid_argument(
            "Observer state contains non-finite values");
    }

    if (!std::isfinite(options_.toleranceSeconds)
        || !(options_.toleranceSeconds > 0.0) ||
        options_.maxIterations == 0) {
        throw std::invalid_argument(
            "Invalid light-time solver options");
    }

    // Initial geometric target position S(t).
    const auto initialTargetState =
        targetProvider.stateAt(receptionEpoch);

    if (!initialTargetState.allFinite()) {
        throw std::runtime_error(
            "Target provider returned a non-finite state");
    }

    Eigen::Vector3d initialLineOfSight =
        initialTargetState.positionKm -
        observerStateAtReception.positionKm;

    double lightTime =
        initialLineOfSight.stableNorm() /
        speedOfLightKmPerSecond;

    if (!(lightTime > 0.0) || !std::isfinite(lightTime)) {
        throw std::runtime_error(
            "Initial light time is not finite");
    }

    for (std::size_t iteration = 1;
         iteration <= options_.maxIterations;
         ++iteration) {

        const TdbEpoch emissionEpoch{
            receptionEpoch.secondsPastJ2000 - lightTime
        };
        if (!std::isfinite(emissionEpoch.secondsPastJ2000)) {
            throw std::runtime_error("Emission epoch is not finite");
        }

        const auto targetState =
            targetProvider.stateAt(emissionEpoch);

        if (!targetState.allFinite()) {
            throw std::runtime_error(
                "Target provider returned a non-finite state");
        }

        // T(t) = S(t - tau) - R(t)
        const Eigen::Vector3d lineOfSight =
            targetState.positionKm -
            observerStateAtReception.positionKm;

        const double rangeKm = lineOfSight.stableNorm();

        if (!(rangeKm > 0.0) || !std::isfinite(rangeKm)) {
            throw std::runtime_error(
                "Invalid observer-target range");
        }

        const double updatedLightTime =
            rangeKm / speedOfLightKmPerSecond;

        const double residual =
            std::abs(updatedLightTime - lightTime);

        if (residual <= options_.toleranceSeconds) {
            return LightTimeSolution{
                .receptionEpoch = receptionEpoch,
                .emissionEpoch = emissionEpoch,
                .lightTimeSeconds = lightTime,
                .residualSeconds = residual,
                .targetStateAtEmission = targetState,
                .geometricLineOfSightKm = lineOfSight,
                .iterations = iteration
            };
        }

        lightTime = updatedLightTime;
    }

    throw std::runtime_error(
        "Light-time iteration did not converge");
}

} // namespace fd::opnav
