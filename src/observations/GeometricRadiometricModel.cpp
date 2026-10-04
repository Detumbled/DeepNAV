#include "observations/GeometricRadiometricModel.hpp"
#include <cmath>
#include <stdexcept>

namespace fd::observations {

fd::filters::MeasurementPrediction
geometricRadiometricPrediction(const Eigen::VectorXd& state,
                               const fd::dynamics::CartesianState& station, LinkDirection direction,
                               double groundBiasSeconds, double groundFractionalFrequency) {
    constexpr double c = 299792.458; // km/s
    if ((state.size() != 6 && state.size() != 8) || !state.allFinite() || !station.allFinite() ||
        !std::isfinite(groundBiasSeconds) || !std::isfinite(groundFractionalFrequency) ||
        (direction != LinkDirection::Uplink && direction != LinkDirection::Downlink))
        throw std::invalid_argument("Invalid geometric radiometric state, station or clock.");
    const Eigen::Vector3d relative = state.head<3>() - station.positionKm;
    const Eigen::Vector3d velocity = state.segment<3>(3) - station.velocityKmPerSec;
    const double range = relative.norm();
    if (!(range > 0.0) || !std::isfinite(range))
        throw std::invalid_argument("Geometric radiometric range must be finite and nonzero.");
    const Eigen::Vector3d unit = relative / range;
    const double rate = unit.dot(velocity);
    fd::filters::MeasurementPrediction result{Eigen::VectorXd(2),
                                              Eigen::MatrixXd::Zero(2, state.size())};
    result.value << range, rate;
    result.jacobian.block<1, 3>(0, 0) = unit.transpose();
    result.jacobian.block<1, 3>(1, 0) = ((velocity - rate * unit) / range).transpose();
    result.jacobian.block<1, 3>(1, 3) = unit.transpose();
    const double sign = direction == LinkDirection::Uplink ? 1.0 : -1.0;
    const double bias = state.size() == 8 ? state[6] : 0.0;
    const double frequency = state.size() == 8 ? state[7] : 0.0;
    result.value[0] += sign * c * (bias - groundBiasSeconds);
    result.value[1] += sign * c * (frequency - groundFractionalFrequency);
    if (state.size() == 8) {
        result.jacobian(0, 6) = sign * c;
        result.jacobian(1, 7) = sign * c;
    }
    return result;
}

} // namespace fd::observations
