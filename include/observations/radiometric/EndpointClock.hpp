#pragma once
#include "observations/radiometric/TwoWayLink.hpp"

namespace fd::observations::radiometric {
// Reading relative to referenceTdb: (t-referenceTdb) + bias(t) - offset(t).
// bias is clock minus its nominal reference; offset is TDB minus that reference.
// An empty offset denotes an ideal TDB reference, not spacecraft proper time.
// Callbacks must be deterministic, continuous and define an increasing clock.
struct EndpointClock {
    double referenceTdb{};
    std::function<double(double)> biasSeconds;
    ClockOffset coordinateMinusReference;
    [[nodiscard]] double correction(double tdb) const;
    [[nodiscard]] long double reading(double tdb) const;
    [[nodiscard]] double coordinateEpoch(long double readingSeconds) const;
};
} // namespace fd::observations::radiometric
