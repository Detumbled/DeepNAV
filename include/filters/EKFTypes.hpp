#pragma once

#include <Eigen/Core>
#include <functional>

namespace fd::filters {

// Cartesian km, km/s; optional clock bias in seconds, fractional frequency and dimensionless SRP
// scale.
enum class StateLayout { Orbit = 6, OrbitClock = 8, OrbitClockSrp = 9 };

struct StatePrediction {
    Eigen::VectorXd state;
    Eigen::MatrixXd transition;
    Eigen::MatrixXd processCovariance;
};

struct MeasurementPrediction {
    Eigen::VectorXd value;
    Eigen::MatrixXd jacobian;
};

// Models are evaluated at the current estimate, never at the simulated truth.
using PropagationFunction =
    std::function<StatePrediction(double fromEpoch, double toEpoch, const Eigen::VectorXd& state)>;
using MeasurementFunction =
    std::function<MeasurementPrediction(double epoch, const Eigen::VectorXd& state)>;

struct InnovationDiagnostics {
    Eigen::VectorXd innovation;
    Eigen::MatrixXd covariance;
    Eigen::VectorXd whitenedInnovation;
    double normalizedInnovationSquared{0.0};
    bool accepted{false};
};

} // namespace fd::filters
