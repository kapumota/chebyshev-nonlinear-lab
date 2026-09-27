#pragma once
#include "../core/matrix.hpp"
#include "../core/vector.hpp"

#ifdef NUMA_HAS_OPENMP
#include <omp.h>
#endif

namespace numa {

// Jacobi paralelo. Reemplaza a gauss_seidel + tri_sup cuando
// se requiere paralelismo. Converge si A es diagonal dominante.
template <std::size_t N>
Vec<N> jacobi_omp(const Matrix<N, N>& A, const Vec<N>& b,
                  std::size_t max_iter = 500,
                  double tol = 1e-12)
{
    Vec<N> x{};
    Vec<N> x_new{};

    for (std::size_t it = 0; it < max_iter; ++it) {
        double diff = 0.0;

        #ifdef NUMA_HAS_OPENMP
        #pragma omp parallel for reduction(+:diff)
        #endif
        for (std::size_t i = 0; i < N; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j)
                if (j != i) s += A[i][j] * x[j];
            x_new[i] = (b[i] - s) / A[i][i];
            diff += (x_new[i] - x[i]) * (x_new[i] - x[i]);
        }

        x = x_new;
        if (std::sqrt(diff) < tol) break;
    }
    return x;
}

} // namespace numa
