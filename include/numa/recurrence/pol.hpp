#pragma once
#include <cmath>

namespace numa {

struct PolResult { double t; double et; };

inline PolResult pol(double t, double L, double a) {
    const double previous_t = t;
    const double p = 0.5 * L * t * t - t + a;
    const double dp = L * t - 1.0;
    t = t - p / dp;
    return { t, std::abs(t - previous_t) };
}

} // namespace numa
