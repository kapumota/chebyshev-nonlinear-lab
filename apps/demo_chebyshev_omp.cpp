#include <numa/methods/chebyshev_omp.hpp>
#include <numa/problems/chandrasekhar.hpp>
#include <cstdio>

#ifdef NUMA_HAS_OPENMP
#include <omp.h>
#endif

int main() {
    using namespace numa;
    using namespace numa::problems;

#ifdef NUMA_HAS_OPENMP
    std::printf("OpenMP habilitado. Hilos disponibles: %d\n", omp_get_max_threads());
#else
    std::printf("OpenMP NO habilitado (modo secuencial)\n");
#endif

    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;

    auto r = chebyshev_omp<8>(f, x0, 1e-12, 20);

    std::printf("CHEBYSHEV OMP (8x8)\n");
    for (std::size_t i = 0; i < 8; ++i)
        std::printf("x[%zu] = %.17g\n", i, r.x[i]);
    std::printf("iteraciones = %zu, error = %.3e\n", r.iterations, r.er_abs);
    return 0;
}
