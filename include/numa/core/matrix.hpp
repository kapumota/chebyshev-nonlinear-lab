#pragma once
#include "vector.hpp"

namespace numa {

template <std::size_t M, std::size_t N>
using Matrix = std::array<std::array<double, N>, M>;

template <std::size_t M, std::size_t N>
Vec<M> matvec(const Matrix<M, N>& A, const Vec<N>& x) {
    Vec<M> y{};
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < N; ++j)
            y[i] += A[i][j] * x[j];
    return y;
}

template <std::size_t N>
Matrix<N, N> identity() {
    Matrix<N, N> I{};
    for (std::size_t i = 0; i < N; ++i) I[i][i] = 1.0;
    return I;
}

template <std::size_t M, std::size_t N>
double norm_inf(const Matrix<M, N>& A) {
    double m = 0.0;
    for (const auto& row : A) {
        double s = 0.0;
        for (auto x : row) s += std::abs(x);
        m = std::max(m, s);
    }
    return m;
}

} // namespace numa
