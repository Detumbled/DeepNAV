#pragma once

#include "filters/EKFTypes.hpp"
#include "filters/filter.hpp"
#include <limits>

namespace fd::filters {

class EKF final : public Filter {
public:
    explicit EKF(StateLayout layout = StateLayout::Orbit);

    void setInitialState(const Eigen::VectorXd& state, const Eigen::MatrixXd& covariance,
                         double epoch) override;
    // Forward prediction commits only after all model outputs have been checked.
    void predictTo(double epoch, const PropagationFunction& propagate);
    // Measurements must be referenced to the current epoch. Latency is a driver concern.
    [[nodiscard]] InnovationDiagnostics
    update(double epoch, const Eigen::VectorXd& observed, const Eigen::MatrixXd& covariance,
           const MeasurementFunction& model,
           double nisLimit = std::numeric_limits<double>::infinity());

    // Joint update with a frozen linearization, including correlated measurement noise.
    void processBatch(const Eigen::VectorXd& residuals, const Eigen::MatrixXd& designMatrix,
                      const Eigen::MatrixXd& measurementCovariance) override;

    [[nodiscard]] StateLayout layout() const noexcept { return layout_; }
    [[nodiscard]] const InnovationDiagnostics& lastInnovation() const noexcept { return last_; }

private:
    [[nodiscard]] InnovationDiagnostics correct(const Eigen::VectorXd& residuals,
                                                const Eigen::MatrixXd& jacobian,
                                                const Eigen::MatrixXd& covariance, double nisLimit);
    StateLayout layout_;
    InnovationDiagnostics last_;
};

} // namespace fd::filters
