#pragma once
#include "newton.hpp"
#include "../linalg/jacobi_omp.hpp"

#ifdef NUMA_HAS_OPENMP
#include <omp.h>
#endif

namespace numa {

// Chebyshev paralelo con OpenMP.
// - eval/jacobian/hessian se paralelizan dentro de la implementación de IFunction.
// - El cálculo de w[i] se paraleliza sobre i.
// - El paso lineal usa Jacobi paralelo.
template <std::size_t N>
SolverResult<N> chebyshev_omp(const IFunction<N>& f,
                              Vec<N> x,
                              double tol = 1e-12,
                              std::size_t max_iter = 20,
                              std::size_t jacobi_max = 500)
{
    double er = 0.0;
    std::size_t k = 0;

    while (k < max_iter) {
        const Vec<N> previous_x = x;
        auto F = f.eval(x);
        auto J = f.jacobian(x);
        auto H = f.hessian(x);

        Vec<N> z = jacobi_omp(J, F, jacobi_max, tol * 1e-2);

        Vec<N> w{};
        #ifdef NUMA_HAS_OPENMP
        #pragma omp parallel for
        #endif
        for (std::size_t i = 0; i < N; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j)
                for (std::size_t kk = 0; kk < N; ++kk)
                    s += H[i][j * N + kk] * z[j] * z[kk];
            w[i] = s;
        }

        Vec<N> u = jacobi_omp(J, w, jacobi_max, tol * 1e-2);

        Vec<N> dx{};
        #ifdef NUMA_HAS_OPENMP
        #pragma omp parallel for
        #endif
        for (std::size_t i = 0; i < N; ++i)
            dx[i] = z[i] + 0.5 * u[i];

        x = x - dx;
        er = norm_inf(x - previous_x);
        ++k;

        if (er < tol) break;
    }
    return { x, k, er, 0.0 };
}

} // namespace numa
