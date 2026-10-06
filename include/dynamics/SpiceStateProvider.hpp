#pragma once
#include "dynamics/StateProvider.hpp"
#include <string>

namespace fd::dynamics {
// Uses the caller-owned CSPICE kernel pool. Geometric J2000/SSB states only.
class SpiceStateProvider final : public StateProvider {
  public:
    explicit SpiceStateProvider(std::string target);
    [[nodiscard]] CartesianState stateAt(opnav::TdbEpoch epoch) const override;
    [[nodiscard]] std::string_view referenceFrame() const noexcept override { return "J2000"; }
    [[nodiscard]] std::string_view referenceOrigin() const noexcept override { return "SSB"; }
    [[nodiscard]] const std::string& target() const noexcept { return target_; }
  private:
    std::string target_;
};
} // namespace fd::dynamics
