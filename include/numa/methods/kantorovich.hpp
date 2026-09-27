#pragma once
#include "newton.hpp"
#include "../recurrence/pol.hpp"

namespace numa {

// Algoritmo 4: Newton-Kantorovich con doble criterio de parada.
// err_abs = ||x_{n+1} - x_n||_inf
// err_teo = |t_{n+1} - t_n|  (sucesión del polinomio p(t))
template <std::size_t N>
SolverResult<N> kantorovich(const IFunction<N>& f,
                            Vec<N> x,
                            double L, double a,
                            double tol = 1e-12,
                            std::size_t max_iter = 20)
{
    const Vec<N> x_ref = x;
    double t = 0.0;
    double et = 1.0;
    double er = 0.0;
    std::size_t k = 0;

    while (et >= tol && k <= max_iter) {
        auto F = f.eval(x);
        auto J = f.jacobian(x);
        gauss_seidel(J, F);
        Vec<N> dx = tri_sup(J, F);

        auto pr = pol(t, L, a);
        t = pr.t; et = pr.et;

        x = x - dx;
        er = norm_inf(x - x_ref);
        ++k;
    }
    return { x, k, er, et };
}

} // namespace numa
