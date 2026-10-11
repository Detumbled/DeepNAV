#include "observations/radiometric/VlbiMeasurement.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>
namespace fd::observations::radiometric {
VlbiMeasurement::VlbiMeasurement(VlbiLink link, const dynamics::StateProvider& one,
    const dynamics::StateProvider& two, double reception, TrajectoryFactory trajectory,
    Eigen::VectorXd steps):link_(std::move(link)),one_(one),two_(two),reception_(reception),
    trajectory_(std::move(trajectory)),steps_(std::move(steps)) {
    if (!std::isfinite(reception) || !trajectory_ || steps_.size()==0 ||
        !steps_.allFinite() || (steps_.array()<=0).any())
        throw std::invalid_argument("Invalid VLBI measurement configuration");
}
double VlbiMeasurement::value(double epoch, const Eigen::VectorXd& state) const {
    if (!std::isfinite(epoch) || state.size()!=steps_.size() || !state.allFinite())
        throw std::invalid_argument("Invalid VLBI measurement state");
    const auto provider=trajectory_(epoch,state);
    if (!provider) throw std::runtime_error("VLBI trajectory factory returned no provider");
    return link_.solve(reception_,*provider,one_,two_).differentialRangeKm;
}
filters::MeasurementPrediction VlbiMeasurement::operator()(double epoch,
    const Eigen::VectorXd& state) const {
    filters::MeasurementPrediction result{Eigen::VectorXd::Constant(1,value(epoch,state)),
                                          Eigen::MatrixXd(1,state.size())};
    for (Eigen::Index j=0;j<state.size();++j) {
        Eigen::VectorXd a=state,b=state;a[j]+=steps_[j];b[j]-=steps_[j];
        result.jacobian(0,j)=(value(epoch,a)-value(epoch,b))/(2*steps_[j]);
    }
    return result;
}
} // namespace fd::observations::radiometric
