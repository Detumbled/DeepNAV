#include "observations/radiometric/OneWayMeasurement.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>
namespace fd::observations::radiometric {
OneWayMeasurement::OneWayMeasurement(OneWayLink link, const dynamics::StateProvider &station,
                                     long double tag, double count, double anchor,
                                     FrequencyRamp ramp, double referenceHz, EndpointClock ground,
                                     ClockOffset nominal, double drift,
                                     TrajectoryFactory trajectory, Eigen::VectorXd steps)
    : link_(std::move(link)), station_(station), tag_(tag), count_(count), clockReference_(anchor),
      referenceHz_(referenceHz), drift_(drift), ramp_(ramp), ground_(std::move(ground)),
      nominal_(std::move(nominal)), trajectory_(std::move(trajectory)), steps_(std::move(steps)) {
    if (!trajectory_ || (steps_.size() != 8 && steps_.size() != 9) || !steps_.allFinite() ||
        (steps_.array() <= 0).any() || !std::isfinite(drift_) || !std::isfinite(tag_) ||
        !std::isfinite(anchor) || !std::isfinite(count) || count <= 0)
        throw std::invalid_argument("Invalid one-way clock measurement configuration.");
}
Eigen::Vector2d OneWayMeasurement::value(double epoch, const Eigen::VectorXd &x) const {
    if (!std::isfinite(epoch) || x.size() != steps_.size() || !x.allFinite())
        throw std::invalid_argument("One-way measurement requires finite orbit-clock state.");
    const auto provider = trajectory_(epoch, x);
    if (!provider)
        throw std::runtime_error("Trajectory factory returned no provider.");
    const EndpointClock onboard{clockReference_,
                                [=, this](double t) {
                                    const double dt = t - epoch;
                                    return x[6] + x[7] * dt + .5 * drift_ * dt * dt;
                                },
                                nominal_};
    const auto count =
        link_.count(tag_, count_, station_, *provider, ground_, onboard, ramp_, referenceHz_);
    return {count.pseudorangeKm, count.residualHz};
}
filters::MeasurementPrediction OneWayMeasurement::operator()(double epoch,
                                                             const Eigen::VectorXd &x) const {
    filters::MeasurementPrediction p{value(epoch, x), Eigen::MatrixXd(2, x.size())};
    for (Eigen::Index j = 0; j < x.size(); ++j) {
        Eigen::VectorXd plus = x, minus = x;
        plus[j] += steps_[j];
        minus[j] -= steps_[j];
        p.jacobian.col(j) = (value(epoch, plus) - value(epoch, minus)) / (2 * steps_[j]);
    }
    return p;
}
} // namespace fd::observations::radiometric
