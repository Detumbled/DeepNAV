#pragma once

#include "Clocks/ClockModel.hpp"
#include "RKF45Integrator.hpp"
#include "filters/EKFTypes.hpp"
#include "perturbations/ForceModel.hpp"

namespace fd::filters {

using fd::perturbations::AccelerationEvaluation;
using fd::perturbations::AccelerationFunction;

struct CartesianPropagationConfig {
    od::RKF45Integrator::Options integrator;
    // Numerical scales for state and STM; entries have Cartesian km and km/s units.
    Eigen::Matrix<double, 6, 1> stateScales{
        (Eigen::Matrix<double, 6, 1>() << 10000, 10000, 10000, 10, 10, 10).finished()};
    // Inertial white acceleration diffusion, km^2/s^3.
    Eigen::Matrix3d accelerationDiffusion{Eigen::Matrix3d::Zero()};
    fd::clocks::ClockParameters clock;
};

// Existing RKF45 integrates scaled orbit/STM/noise equations; clock dynamics are exact.
class CartesianPropagator {
public:
    explicit CartesianPropagator(AccelerationFunction acceleration,
                                 CartesianPropagationConfig config = {});
    [[nodiscard]] StatePrediction operator()(double fromEpoch, double toEpoch,
                                             const Eigen::VectorXd& state) const;

private:
    AccelerationFunction acceleration_;
    CartesianPropagationConfig config_;
    fd::clocks::ClockModel clock_;
    od::RKF45Integrator integrator_;
};

} // namespace fd::filters
