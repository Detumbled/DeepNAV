#pragma once
#include "RKF45Integrator.hpp"
#include "dynamics/StateProvider.hpp"
#include "perturbations/ForceModel.hpp"
#include <string>

namespace fd::dynamics {
// Numerical trajectory with Hermite dense output; no CSPICE dependency.
// Acceleration receives absolute TDB epochs. An optional moving-origin provider
// converts relative states to its fixed origin; it must outlive this object.
class IntegratedStateProvider final : public StateProvider {
  public:
    IntegratedStateProvider(double initialTdb, const Eigen::Matrix<double, 6, 1> &initial,
                            double endTdb, const perturbations::AccelerationFunction &acceleration,
                            od::RKF45Integrator::Options options = {}, std::string frame = "J2000",
                            std::string origin = "SSB",
                            const StateProvider *movingOrigin = nullptr);
    [[nodiscard]] CartesianState stateAt(opnav::TdbEpoch epoch) const override;
    [[nodiscard]] CartesianState relativeStateAt(double epochTdb) const;
    [[nodiscard]] std::string_view referenceFrame() const noexcept override { return frame_; }
    [[nodiscard]] std::string_view referenceOrigin() const noexcept override { return origin_; }

  private:
    double initialTdb_;
    od::EphemerisInterpolator ephemeris_;
    std::string frame_, origin_;
    const StateProvider *movingOrigin_;
};
} // namespace fd::dynamics
