#include <cstdio>
#include <cuda_runtime.h>
#include <numa/cuda/chebyshev_kernel.cuh>
#include <numa/functions/gauss_legendre.hpp>

constexpr std::size_t N = 8;

int main() {
    using namespace numa;
    using namespace numa::cuda;

    // Inicializar A con Gauss-Legendre
    ChandrasekharDevice<N> f;
    f.init_A(GaussLegendre<N>::t().data(), GaussLegendre<N>::w().data());

    // Estado inicial
    double x_host[N];
    for (std::size_t i = 0; i < N; ++i) x_host[i] = 1.0;

    // Lanzar
    double x_out[N];
    int iter = 0;
    double er = 0.0;

    chebyshev_gpu<ChandrasekharDevice<N>, N>(
        f, x_host, x_out, &iter, &er,
        /*tol=*/1e-12, /*max_iter=*/20, /*jacobi_max=*/200);

    std::printf("CHEBYSHEV GPU (N=%zu)\n", N);
    for (std::size_t i = 0; i < N; ++i)
        std::printf("x[%zu] = %.17g\n", i, x_out[i]);
    std::printf("iteraciones = %d, error = %.3e\n", iter, er);

    return 0;
}
