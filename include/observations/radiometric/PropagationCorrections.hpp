#pragma once
#include "observations/radiometric/TwoWayLink.hpp"

namespace fd::observations::radiometric {
// Endpoint/body states must use the same inertial frame and origin.
[[nodiscard]] double pointMassShapiroSeconds(const LegSolution& leg,
                                             const dynamics::StateProvider& body, double gmKm3S2);
// Simple non-dispersive zenith-path mapping. Use only above an elevation mask;
// measured wet/dry meteorology and mapping functions can replace this model.
[[nodiscard]] double zenithPathDelaySeconds(double zenithPathM, double elevationRadians);
} // namespace fd::observations::radiometric
