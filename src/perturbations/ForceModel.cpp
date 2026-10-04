#include "perturbations/ForceModel.hpp"
#include <stdexcept>
#include <utility>

namespace fd::perturbations {

AccelerationFunction sumAccelerations(std::vector<AccelerationFunction> forces) {
    for (const auto& force : forces)
        if (!force)
            throw std::invalid_argument("Force composition contains an empty model.");
    return [forces = std::move(forces)](double epoch, const Eigen::Vector3d& position,
                                        const Eigen::Vector3d& velocity) {
        AccelerationEvaluation result;
        for (const auto& force : forces) {
            const auto contribution = force(epoch, position, velocity);
            result.acceleration += contribution.acceleration;
            result.positionJacobian += contribution.positionJacobian;
            result.velocityJacobian += contribution.velocityJacobian;
        }
        return result;
    };
}

} // namespace fd::perturbations
