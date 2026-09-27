#pragma once
#include "../core/vector.hpp"
#include "../core/matrix.hpp"
#include "../functions/function.hpp"
#include "../linalg/gauss_seidel.hpp"
#include "../linalg/tri_sup.hpp"

namespace numa {

template <std::size_t N>
struct SolverResult {
    Vec<N> x;
    std::size_t iterations;
    double er_abs;
    double er_teorico;
};

// Algoritmo 1 del Anexo: Newton clásico.
// x_{n+1} = x_n - J(x_n)^{-1} F(x_n)
template <std::size_t N>
SolverResult<N> newton(const IFunction<N>& f,
                       Vec<N> x,
                       double tol = 1e-12,
                       std::size_t max_iter = 20)
{
    double er = 0.0;
    std::size_t k = 0;

    while (k < max_iter) {
        const Vec<N> previous_x = x;
        auto J = f.jacobian(x);
        auto F = f.eval(x);
        gauss_seidel(J, F);
        const Vec<N> dx = tri_sup(J, F);
        x = x - dx;
        er = norm_inf(x - previous_x);
        ++k;

        if (er < tol) break;
    }

    return { x, k, er, 0.0 };
}

} // namespace numa
