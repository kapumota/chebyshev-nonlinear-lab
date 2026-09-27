#pragma once
#include <cuda_runtime.h>
#include <cstddef>
#include <cmath>

namespace numa::cuda {

// Jacobi iterativo en device. Un bloque por sistema, N hilos.
// Reemplaza a Gauss-Seidel (secuencial) en GPU.
template <std::size_t N>
__global__
void jacobi_kernel(const double* A, const double* b, double* x,
                   int max_iter, double tol)
{
    __shared__ double xs[N];
    __shared__ double xnew[N];

    const int i = threadIdx.x;
    if (i >= (int)N) return;

    xs[i] = 0.0;
    __syncthreads();

    for (int it = 0; it < max_iter; ++it) {
        double s = 0.0;
        for (std::size_t j = 0; j < N; ++j)
            if ((int)j != i) s += A[i * N + j] * xs[j];
        xnew[i] = (b[i] - s) / A[i * N + i];
        __syncthreads();

        double diff = xnew[i] - xs[i];
        xs[i] = xnew[i];
        __syncthreads();

        if (fabs(diff) < tol) break;
    }

    x[i] = xs[i];
}

// Wrapper host
template <std::size_t N>
void jacobi_gpu(const double* d_A, const double* d_b, double* d_x,
                int max_iter = 500, double tol = 1e-10)
{
    jacobi_kernel<N><<<1, N>>>(d_A, d_b, d_x, max_iter, tol);
    cudaDeviceSynchronize();
}

} // namespace numa::cuda
