#pragma once
#include <cmath>

namespace numa {

struct PolResult { double t; double et; };

inline PolResult pol(double t, double L, double a) {
    const double t1 = t;
    t = t - (L * t * t - t + a) / (L * t - 1.0);
    return { t, std::abs(t - t1) };
}

} // namespace numa
