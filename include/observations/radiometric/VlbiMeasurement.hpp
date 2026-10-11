#pragma once
#include "filters/EKFTypes.hpp"
#include "observations/radiometric/TrajectoryFactory.hpp"
#include "observations/radiometric/VlbiLink.hpp"

namespace fd::observations::radiometric {
class VlbiMeasurement {
  public:
    VlbiMeasurement(VlbiLink link, const dynamics::StateProvider& stationOne,
        const dynamics::StateProvider& stationTwo, double receptionTdb,
        TrajectoryFactory trajectory, Eigen::VectorXd steps);
    [[nodiscard]] double value(double stateEpochTdb, const Eigen::VectorXd& state) const;
    [[nodiscard]] filters::MeasurementPrediction operator()(double stateEpochTdb,
        const Eigen::VectorXd& state) const;
  private:
    VlbiLink link_;
    const dynamics::StateProvider& one_;
    const dynamics::StateProvider& two_;
    double reception_;
    TrajectoryFactory trajectory_;
    Eigen::VectorXd steps_;
};
} // namespace fd::observations::radiometric
