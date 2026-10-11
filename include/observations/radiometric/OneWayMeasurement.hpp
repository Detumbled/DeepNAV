#pragma once
#include "Clocks/Types.hpp"
#include "filters/EKFTypes.hpp"
#include "observations/radiometric/OneWayLink.hpp"
#include "observations/radiometric/TrajectoryFactory.hpp"

namespace fd::observations::radiometric {
// Raw spacecraft time tag relative to receiverReferenceTdb; never pre-corrected
// with truth. Clock bias/frequency belong to state[6:8] at stateEpochTdb.
class OneWayMeasurement {
  public:
    OneWayMeasurement(OneWayLink link, const dynamics::StateProvider &station,
                      long double receptionClockSeconds, double countClockSeconds,
                      double receiverReferenceTdb, FrequencyRamp ramp, double referenceHz,
                      EndpointClock ground, ClockOffset receiverNominalOffset,
                      double frequencyDriftPerSecond, TrajectoryFactory trajectory,
                      Eigen::VectorXd finiteDifferenceSteps);
    [[nodiscard]] Eigen::Vector2d value(double stateEpochTdb, const Eigen::VectorXd &state) const;
    [[nodiscard]] filters::MeasurementPrediction operator()(double stateEpochTdb,
                                                            const Eigen::VectorXd &state) const;

  private:
    OneWayLink link_;
    const dynamics::StateProvider &station_;
    long double tag_;
    double count_, clockReference_, referenceHz_, drift_;
    FrequencyRamp ramp_;
    EndpointClock ground_;
    ClockOffset nominal_;
    TrajectoryFactory trajectory_;
    Eigen::VectorXd steps_;
};
} // namespace fd::observations::radiometric
