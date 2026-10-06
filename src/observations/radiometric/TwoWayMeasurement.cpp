#include "observations/radiometric/TwoWayMeasurement.hpp"
#include <stdexcept>
#include <utility>
namespace fd::observations::radiometric {
TwoWayMeasurement::TwoWayMeasurement(TwoWayLink link, const dynamics::StateProvider &station,
                                     double reception, double count, FrequencyRamp ramp,
                                     double referenceHz, ClockOffset clock,
                                     TrajectoryFactory trajectory, Eigen::VectorXd steps)
    : link_(std::move(link)), station_(station), reception_(reception), count_(count),
      referenceHz_(referenceHz), ramp_(ramp), clock_(std::move(clock)),
      trajectory_(std::move(trajectory)), steps_(std::move(steps)) {
    if (!trajectory_ || steps_.size() == 0 || !steps_.allFinite() || (steps_.array() <= 0).any())
        throw std::invalid_argument(
            "Trajectory factory and positive differentiation steps required.");
}
Eigen::Vector2d TwoWayMeasurement::value(double epoch, const Eigen::VectorXd &state) const {
    if (state.size() != steps_.size() || !state.allFinite())
        throw std::invalid_argument("Invalid state for two-way measurement.");
    const auto provider = trajectory_(epoch, state);
    if (!provider)
        throw std::runtime_error("Trajectory factory returned no provider.");
    const auto count =
        link_.count(reception_, count_, station_, *provider, ramp_, referenceHz_, clock_);
    return {count.end.rangeKm(), count.residualHz};
}
filters::MeasurementPrediction TwoWayMeasurement::operator()(double epoch,
                                                             const Eigen::VectorXd &state) const {
    filters::MeasurementPrediction result{value(epoch, state), Eigen::MatrixXd(2, state.size())};
    for (Eigen::Index j = 0; j < state.size(); ++j) {
        Eigen::VectorXd plus = state, minus = state;
        plus[j] += steps_[j];
        minus[j] -= steps_[j];
        result.jacobian.col(j) = (value(epoch, plus) - value(epoch, minus)) / (2 * steps_[j]);
    }
    return result;
}
} // namespace fd::observations::radiometric
