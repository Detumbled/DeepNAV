#pragma once

#include "dynamics/StateProvider.hpp"

#include <string>
#include <string_view>

namespace fd::dynamics {

class LinearStateProvider final : public StateProvider {
public:
    LinearStateProvider(
        opnav::TdbEpoch referenceEpoch,
        CartesianState referenceState,
        std::string referenceFrame,
        std::string referenceOrigin);

    [[nodiscard]]
    CartesianState stateAt(
        opnav::TdbEpoch epoch) const override;

    [[nodiscard]]
    std::string_view referenceFrame() const noexcept override;

    [[nodiscard]]
    std::string_view referenceOrigin() const noexcept override;

private:
    opnav::TdbEpoch referenceEpoch_;
    CartesianState referenceState_;
    std::string referenceFrame_;
    std::string referenceOrigin_;
};

} // namespace fd::dynamics
