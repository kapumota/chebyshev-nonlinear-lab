#pragma once
#include <cmath>

namespace numa {

inline double sigma(double a) {
    const double s1 = 2.0 * (2.0 * a * a - 3.0 * a - 1.0);
    const double s2 = 2.0 * std::sqrt(1.0 + 8.0 * a - 4.0 * a * a);
    return (s1 + s2) / (a * (1.0 - 2.0 * a));
}

} // namespace numa
