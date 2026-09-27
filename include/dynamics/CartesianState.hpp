#pragma once

#include <Eigen/Core>

namespace fd::dynamics {

struct CartesianState {
    Eigen::Vector3d positionKm {Eigen::Vector3d::Zero()};
    Eigen::Vector3d velocityKmPerSec {Eigen::Vector3d::Zero()};

    [[nodiscard]]
    bool allFinite() const noexcept {
        return positionKm.allFinite()
            && velocityKmPerSec.allFinite();
    }
};

} // namespace fd::dynamics
