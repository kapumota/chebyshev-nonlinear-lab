#pragma once
#include <cuda_runtime.h>
#include <cstddef>
#include <cmath>
#include "chandrasekhar_device.hpp"

namespace numa::cuda {

// Chebyshev completo en GPU. Un bloque por sistema, N hilos.
// Cada iteración:
//   1. F = f(x)
//   2. J = f.jacobian(x)
//   3. H = f.hessian(x)
//   4. z = J^{-1} F       (Jacobi iterativo)
//   5. w[i] = H[i](z,z)
//   6. u = J^{-1} w       (Jacobi iterativo)
//   7. dx = z + 0.5 u
//   8. x -= dx
template <typename F, std::size_t N>
__global__
void chebyshev_kernel(const F* f,
                      const double* x_in,
                      double* x_out,
                      int* iter_out,
                      double* er_out,
                      double tol,
                      int max_iter,
                      int jacobi_max)
{
    __shared__ double x[N];
    __shared__ double F[N];
    __shared__ double J[N * N];
    __shared__ double H[N * N * N];
    __shared__ double z[N];
    __shared__ double w[N];
    __shared__ double u[N];
    __shared__ double dx[N];

    const int i = threadIdx.x;
    if (i >= (int)N) return;

    for (std::size_t k = 0; k < N; ++k) x[k] = x_in[k];
    __syncthreads();

    int it = 0;
    double err = 0.0;

    while (it < max_iter) {
        // 1. F = f(x)
        f->eval(x, F);
        __syncthreads();

        // 2. J = f.jacobian(x)
        f->jacobian(x, J);
        __syncthreads();

        // 3. H = f.hessian(x)
        f->hessian(x, H);
        __syncthreads();

        // 4. z = J^{-1} F  (Jacobi)
        z[i] = 0.0;
        __syncthreads();
        for (int js = 0; js < jacobi_max; ++js) {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j)
                if ((int)j != i) s += J[i * N + j] * z[j];
            double zi_new = (F[i] - s) / J[i * N + i];
            __syncthreads();
            z[i] = zi_new;
            __syncthreads();
        }

        // 5. w[i] = Σ_{j,k} H[i][j*N+k] z[j] z[k]
        {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j)
                for (std::size_t k = 0; k < N; ++k)
                    s += H[i * N * N + j * N + k] * z[j] * z[k];
            w[i] = s;
        }
        __syncthreads();

        // 6. u = J^{-1} w  (Jacobi)
        u[i] = 0.0;
        __syncthreads();
        for (int js = 0; js < jacobi_max; ++js) {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j)
                if ((int)j != i) s += J[i * N + j] * u[j];
            double ui_new = (w[i] - s) / J[i * N + i];
            __syncthreads();
            u[i] = ui_new;
            __syncthreads();
        }

        // 7. dx = z + 0.5 u
        dx[i] = z[i] + 0.5 * u[i];
        __syncthreads();

        // 8. x -= dx, medir error
        double local_err = fabs(dx[i]);
        x[i] -= dx[i];
        __syncthreads();

        // Reducción del error (solo hilo 0)
        if (i == 0) {
            err = 0.0;
            for (std::size_t k = 0; k < N; ++k)
                err = fmax(err, fabs(dx[k]));
        }
        __syncthreads();

        ++it;
        if (err < tol) break;
    }

    x_out[i] = x[i];
    if (i == 0) { *iter_out = it; *er_out = err; }
}

// Wrapper host
template <typename F, std::size_t N>
void chebyshev_gpu(const F& host_f,
                   const double* x_host,
                   double* x_out_host,
                   int* iter_out,
                   double* er_out,
                   double tol = 1e-12,
                   int max_iter = 20,
                   int jacobi_max = 200)
{
    F* d_f = nullptr;
    double* d_x = nullptr;
    int* d_iter = nullptr;
    double* d_er = nullptr;

    cudaMalloc(&d_f, sizeof(F));
    cudaMalloc(&d_x, N * sizeof(double));
    cudaMalloc(&d_iter, sizeof(int));
    cudaMalloc(&d_er, sizeof(double));

    cudaMemcpy(d_f, &host_f, sizeof(F), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x_host, N * sizeof(double), cudaMemcpyHostToDevice);

    chebyshev_kernel<F, N><<<1, (int)N>>>(
        d_f, d_x, d_x, d_iter, d_er, tol, max_iter, jacobi_max);

    cudaError_t err_cuda = cudaDeviceSynchronize();
    if (err_cuda != cudaSuccess)
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err_cuda));

    cudaMemcpy(x_out_host, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(iter_out, d_iter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(er_out, d_er, sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_f);
    cudaFree(d_x);
    cudaFree(d_iter);
    cudaFree(d_er);
}

} // namespace numa::cuda
