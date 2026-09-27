#pragma once
#include "newton.hpp"
#include <limits>

namespace numa {

// Algoritmo 2: Newton con diferencias finitas hacia adelante.
// h_j = sqrt(eps) * max(|x_j|, 1)
template <std::size_t N>
Matrix<N, N> jacobian_fd(const IFunction<N>& f, const Vec<N>& x) {
    Matrix<N, N> J{};
    Vec<N> F0 = f.eval(x);
    for (std::size_t j = 0; j < N; ++j) {
        const double h = std::sqrt(std::numeric_limits<double>::epsilon())
                         * std::max(std::abs(x[j]), 1.0);
        Vec<N> xj = x;
        xj[j] += h;
        Vec<N> Fj = f.eval(xj);
        for (std::size_t i = 0; i < N; ++i)
            J[i][j] = (Fj[i] - F0[i]) / h;
    }
    return J;
}

template <std::size_t N>
SolverResult<N> newton_df(const IFunction<N>& f,
                          Vec<N> x,
                          double tol = 1e-12,
                          std::size_t max_iter = 20)
{
    const Vec<N> x_ref = x;
    double er = 0.0;
    std::size_t k = 0;
    while (k == 0 || (er >= tol && k <= max_iter)) {
        auto J = jacobian_fd(f, x);
        auto F = f.eval(x);
        gauss_seidel(J, F);
        Vec<N> dx = tri_sup(J, F);
        x = x - dx;
        er = norm_inf(x - x_ref);
        ++k;
    }
    return { x, k, er, 0.0 };
}

} // namespace numa
