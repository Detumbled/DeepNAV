#pragma once
#include "filters/EKFTypes.hpp"
#include "observations/radiometric/TwoWayLink.hpp"
#include <memory>
namespace fd::observations::radiometric {
// Build a trajectory around an estimated state at its reference epoch. The
// factory owns propagation/latency choices and must cover all retarded epochs.
using TrajectoryFactory = std::function<std::unique_ptr<dynamics::StateProvider>(
    double stateEpochTdb, const Eigen::VectorXd &state)>;
class TwoWayMeasurement {
  public:
    TwoWayMeasurement(TwoWayLink link, const dynamics::StateProvider &station, double receptionTdb,
                      double countSeconds, FrequencyRamp ramp, double referenceHz,
                      ClockOffset clock, TrajectoryFactory trajectory,
                      Eigen::VectorXd finiteDifferenceSteps);
    [[nodiscard]] Eigen::Vector2d value(double stateEpochTdb, const Eigen::VectorXd &state) const;
    [[nodiscard]] filters::MeasurementPrediction operator()(double stateEpochTdb,
                                                            const Eigen::VectorXd &state) const;

  private:
    TwoWayLink link_;
    const dynamics::StateProvider &station_;
    double reception_, count_, referenceHz_;
    FrequencyRamp ramp_;
    ClockOffset clock_;
    TrajectoryFactory trajectory_;
    Eigen::VectorXd steps_;
};
} // namespace fd::observations::radiometric
