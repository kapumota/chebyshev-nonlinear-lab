// Verificación de consistencia numérica entre uno y varios hilos.
// No mide aceleración y no produce evidencia de rendimiento.

#include <numa/linalg/poisson_mg.hpp>
#include <numa/methods/krylov_newton.hpp>
#include <numa/problems/bratu3d.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>

#ifdef NUMA_HAS_OPENMP
#include <omp.h>
#endif

using namespace numa::op;

namespace {

double relative_error(const DVec& a, const DVec& b) {
    DVec difference = a;
    axpy(-1.0, b, difference);
    const double denominator = norm2(b);
    return denominator > 0.0
        ? norm2(difference) / denominator
        : norm2(difference);
}

} // namespace

int main() {
#ifndef NUMA_HAS_OPENMP
    std::printf("OpenMP no está habilitado; verificación omitida.\n");
    return 0;
#else
    Bratu3D problem(17, 1.0);
    const DVec x0(problem.size(), 0.0);

    KrylovParams params;
    params.recycle = true;

    const int previous = omp_get_max_threads();

    omp_set_num_threads(1);
    const auto one =
        chebyshev_inexact(
            problem, x0, params, PoissonMGProvider(17));

    const int requested = std::min(4, omp_get_num_procs());
    omp_set_num_threads(std::max(1, requested));
    const auto many =
        chebyshev_inexact(
            problem, x0, params, PoissonMGProvider(17));

    omp_set_num_threads(previous);

    const double error = relative_error(many.x, one.x);
    const bool ok =
        one.converged &&
        many.converged &&
        one.iters == many.iters &&
        one.n_F == many.n_F &&
        one.n_JVP == many.n_JVP &&
        one.n_D2 == many.n_D2 &&
        error < 1e-9;

    std::printf(
        "threads=%d iter_1=%zu iter_many=%zu "
        "nJVP_1=%zu nJVP_many=%zu rel=%.3e\n",
        std::max(1, requested),
        one.iters, many.iters,
        one.n_JVP, many.n_JVP, error);

    return ok ? 0 : 1;
#endif
}
