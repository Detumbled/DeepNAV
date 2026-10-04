#pragma once

#include "dynamics/CartesianState.hpp"
#include "filters/EKFTypes.hpp"

namespace fd::observations {

enum class LinkDirection { Uplink, Downlink };

// Simultaneous-epoch range/range-rate benchmark, NOT a light-time/count-time model.
// Station state must use the same inertial frame, origin and epoch as the spacecraft.
[[nodiscard]] fd::filters::MeasurementPrediction geometricRadiometricPrediction(
    const Eigen::VectorXd& state, const fd::dynamics::CartesianState& station,
    LinkDirection direction = LinkDirection::Uplink, double groundBiasSeconds = 0.0,
    double groundFractionalFrequency = 0.0);

} // namespace fd::observations
