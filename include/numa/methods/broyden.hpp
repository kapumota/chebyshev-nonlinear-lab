#pragma once
#include "newton.hpp"
#include <limits>

namespace numa {

// Algoritmo 3: Método de Broyden.
// C_{n+1} = C_n + (y - C_n s) s^T / (s^T s)
template <std::size_t N>
SolverResult<N> broyden(const IFunction<N>& f,
                        Vec<N> x,
                        double tol = 1e-12,
                        std::size_t max_iter = 20)
{
    Matrix<N, N> C = f.jacobian(x);
    double er = 0.0;
    std::size_t k = 0;

    while (k < max_iter) {
        const Vec<N> F = f.eval(x);
        auto D = C;
        auto d = F;
        gauss_seidel(D, d);
        const Vec<N> correction = tri_sup(D, d);
        const Vec<N> step = -1.0 * correction;
        const Vec<N> next_x = x + step;
        const Vec<N> y = f.eval(next_x) - F;

        const Vec<N> Cstep = matvec(C, step);
        const Vec<N> numerator = y - Cstep;
        double denominator = 0.0;
        for (const auto value : step) denominator += value * value;

        if (denominator > std::numeric_limits<double>::epsilon()) {
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < N; ++j)
                    C[i][j] += numerator[i] * step[j] / denominator;
        }

        x = next_x;
        er = norm_inf(step);
        ++k;

        if (er < tol) break;
    }

    return { x, k, er, 0.0 };
}

} // namespace numa
