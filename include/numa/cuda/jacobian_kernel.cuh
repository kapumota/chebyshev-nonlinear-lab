#pragma once
#include <cuda_runtime.h>
#include <cstddef>

namespace numa::cuda {

// Kernel: J[i][j] por diferencias finitas centrales.
// Cada hilo (i,j) calcula una entrada.
// grid = ceil(N/16) x ceil(N/16), block = 16 x 16
template <typename F, std::size_t N>
__global__
void jacobian_fd_kernel(const F* f,
                        const double* x,
                        double* J,
                        double h)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= (int)N || j >= (int)N) return;

    double x_plus[N];
    double x_minus[N];
    double F_plus[N];
    double F_minus[N];

    for (std::size_t k = 0; k < N; ++k) {
        x_plus[k]  = x[k];
        x_minus[k] = x[k];
    }
    x_plus[j]  += h;
    x_minus[j] -= h;

    f->eval(x_plus,  F_plus);
    f->eval(x_minus, F_minus);

    J[i * N + j] = (F_plus[i] - F_minus[i]) / (2.0 * h);
}

// Wrapper host
template <typename F, std::size_t N>
void jacobian_fd_gpu(const F* d_f, const double* d_x, double* d_J, double h)
{
    dim3 block(16, 16);
    dim3 grid(((int)N + 15) / 16, ((int)N + 15) / 16);
    jacobian_fd_kernel<F, N><<<grid, block>>>(d_f, d_x, d_J, h);
    cudaDeviceSynchronize();
}

} // namespace numa::cuda
