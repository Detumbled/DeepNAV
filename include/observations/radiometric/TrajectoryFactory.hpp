#pragma once
#include "dynamics/StateProvider.hpp"
#include <Eigen/Core>
#include <functional>
#include <memory>
namespace fd::observations::radiometric {
using TrajectoryFactory = std::function<std::unique_ptr<dynamics::StateProvider>(
    double stateEpochTdb, const Eigen::VectorXd &state)>;
} // namespace fd::observations::radiometric
