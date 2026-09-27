#pragma once

#include "dynamics/StateProvider.hpp"

#include <string>
#include <string_view>

namespace fd::dynamics {

class ConstantStateProvider final : public StateProvider {
public:
    ConstantStateProvider(
        CartesianState state,
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
    CartesianState state_;
    std::string referenceFrame_;
    std::string referenceOrigin_;
};

} // namespace fd::dynamics
