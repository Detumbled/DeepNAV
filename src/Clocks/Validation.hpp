#pragma once

#include "Clocks/Types.hpp"

#include <cmath>
#include <stdexcept>

namespace fd::clocks::detail {

inline void finite(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("Clock inputs must be finite.");
}

inline void nonnegative(double value) {
    finite(value);
    if (value < 0.0) throw std::invalid_argument("Clock noise and uncertainty must be nonnegative.");
}

inline void positive(double value) {
    finite(value);
    if (value <= 0.0) throw std::invalid_argument("Clock time intervals must be positive.");
}

inline void validate(const ClockParameters& p) {
    finite(p.frequency_drift_per_s);
    nonnegative(p.q_bias_s);
    nonnegative(p.q_frequency_per_s);
}

inline void validate(const ClockState& s) {
    finite(s.bias_s);
    finite(s.fractional_frequency);
}

inline void validate(const ClockCovariance& p) {
    if (!p.allFinite() || p(0, 0) < 0.0 || p(1, 1) < 0.0
        || p(0, 1) != p(1, 0)
        || std::abs(p(0, 1)) > std::sqrt(p(0, 0)) * std::sqrt(p(1, 1)))
        throw std::invalid_argument("Clock covariance must be finite, symmetric and positive semidefinite.");
}

inline double checked(double value) {
    if (!std::isfinite(value)) throw std::overflow_error("Clock calculation overflowed.");
    return value;
}

} // namespace fd::clocks::detail
