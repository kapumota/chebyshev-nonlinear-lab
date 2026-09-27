#pragma once
#include "newton.hpp"

namespace numa {

// Algoritmo 3: Método de Broyden.
// C_{n+1} = C_n + (y - C_n s) s^T / (s^T s)
template <std::size_t N>
SolverResult<N> broyden(const IFunction<N>& f,
                        Vec<N> x,
                        double tol = 1e-12,
                        std::size_t max_iter = 20)
{
    const Vec<N> x_ref = x;
    Matrix<N, N> C = f.jacobian(x);

    auto F0 = f.eval(x);
    auto D = C;
    auto d = F0;
    gauss_seidel(D, d);
    Vec<N> s = tri_sup(D, d);

    double er = norm_inf(s);
    std::size_t k = 0;

    while (er >= tol && k <= max_iter) {
        Vec<N> xx = x - s;
        Vec<N> y  = f.eval(xx) - f.eval(x);

        Vec<N> Cs = matvec(C, s);
        Vec<N> num = y - Cs;
        double den = 0.0;
        for (auto v : s) den += v * v;

        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                C[i][j] += num[i] * s[j] / den;

        x = xx;
        auto F = f.eval(x);
        D = C; d = F;
        gauss_seidel(D, d);
        s = tri_sup(D, d);
        er = norm_inf(x - x_ref);
        ++k;
    }
    return { x, k, er, 0.0 };
}

} // namespace numa
