#pragma once
#include "newton.hpp"
#include "../linalg/jacobi_omp.hpp"

namespace numa {

// Newton paralelo con Jacobi OpenMP.
// Requiere que el Jacobiano sea diagonal dominante para converger.
// Para sistemas densos pequeños (N<=16) usar newton() secuencial.
template <std::size_t N>
SolverResult<N> newton_omp(const IFunction<N>& f,
                           Vec<N> x,
                           double tol = 1e-12,
                           std::size_t max_iter = 20,
                           std::size_t jacobi_max = 500)
{
    double er = 0.0;
    std::size_t k = 0;
    while (k < max_iter) {
        const Vec<N> previous_x = x;
        auto J = f.jacobian(x);
        auto F = f.eval(x);
        Vec<N> dx = jacobi_omp(J, F, jacobi_max, tol * 1e-2);
        x = x - dx;
        er = norm_inf(x - previous_x);
        ++k;

        if (er < tol) break;
    }
    return { x, k, er, 0.0 };
}

} // namespace numa
