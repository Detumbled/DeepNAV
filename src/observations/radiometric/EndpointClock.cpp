#include "observations/radiometric/EndpointClock.hpp"
#include <cmath>
#include <stdexcept>

namespace fd::observations::radiometric {
double EndpointClock::correction(double t) const {
    if (!std::isfinite(t) || !std::isfinite(referenceTdb))
        throw std::invalid_argument("Clock epochs must be finite.");
    const double result = (biasSeconds ? biasSeconds(t) : 0) -
                          (coordinateMinusReference ? coordinateMinusReference(t) : 0);
    if (!std::isfinite(result))
        throw std::runtime_error("Non-finite endpoint clock correction.");
    return result;
}
long double EndpointClock::reading(double t) const {
    return (static_cast<long double>(t) - referenceTdb) + correction(t);
}
double EndpointClock::coordinateEpoch(long double tag) const {
    if (!std::isfinite(tag) || !std::isfinite(referenceTdb))
        throw std::invalid_argument("Clock tag must be finite.");
    long double elapsed = tag;
    for (int i = 0; i < 32; ++i) {
        const double t = static_cast<double>(static_cast<long double>(referenceTdb) + elapsed);
        const long double next = tag - correction(t);
        if (std::abs(next - elapsed) <= 1e-12L)
            return static_cast<double>(static_cast<long double>(referenceTdb) + next);
        elapsed = next;
    }
    throw std::runtime_error("Endpoint clock inversion did not converge.");
}
} // namespace fd::observations::radiometric
