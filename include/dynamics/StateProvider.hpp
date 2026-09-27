#pragma once

#include "dynamics/CartesianState.hpp"
#include "opnav/Types.hpp"

#include <string_view>

namespace fd::dynamics {

class StateProvider {
public:
    virtual ~StateProvider() = default;

    [[nodiscard]]
    virtual CartesianState stateAt(
        opnav::TdbEpoch epoch) const = 0;

    [[nodiscard]]
    virtual std::string_view referenceFrame() const noexcept = 0;

    [[nodiscard]]
    virtual std::string_view referenceOrigin() const noexcept = 0;
};

} // namespace fd::dynamics
