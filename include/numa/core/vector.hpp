#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <algorithm>

namespace numa {

template <std::size_t N>
using Vec = std::array<double, N>;

template <std::size_t N>
double norm_inf(const Vec<N>& v) {
    double m = 0.0;
    for (auto x : v) m = std::max(m, std::abs(x));
    return m;
}

template <std::size_t N>
double norm_2(const Vec<N>& v) {
    double s = 0.0;
    for (auto x : v) s += x * x;
    return std::sqrt(s);
}

template <std::size_t N>
Vec<N> operator-(const Vec<N>& a, const Vec<N>& b) {
    Vec<N> r;
    for (std::size_t i = 0; i < N; ++i) r[i] = a[i] - b[i];
    return r;
}

template <std::size_t N>
Vec<N> operator+(const Vec<N>& a, const Vec<N>& b) {
    Vec<N> r;
    for (std::size_t i = 0; i < N; ++i) r[i] = a[i] + b[i];
    return r;
}

template <std::size_t N>
Vec<N> operator*(double s, const Vec<N>& v) {
    Vec<N> r;
    for (std::size_t i = 0; i < N; ++i) r[i] = s * v[i];
    return r;
}

template <std::size_t N>
Vec<N> operator/(const Vec<N>& v, double s) {
    Vec<N> r;
    for (std::size_t i = 0; i < N; ++i) r[i] = v[i] / s;
    return r;
}

} // namespace numa
