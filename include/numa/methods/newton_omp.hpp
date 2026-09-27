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
    const Vec<N> x_ref = x;
    double er = 0.0;
    std::size_t k = 0;
    while (k == 0 || (er >= tol && k <= max_iter)) {
        auto J = f.jacobian(x);
        auto F = f.eval(x);
        Vec<N> dx = jacobi_omp(J, F, jacobi_max, tol * 1e-2);
        x = x - dx;
        er = norm_inf(x - x_ref);
        ++k;
    }
    return { x, k, er, 0.0 };
}

} // namespace numa
