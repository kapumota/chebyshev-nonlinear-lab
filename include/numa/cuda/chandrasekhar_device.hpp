#pragma once
#include "function_device.hpp"

// Versión device de Chandrasekhar. El arreglo A[N][N] se copia
// al device con cudaMemcpy. N debe ser conocido en compilación.

namespace numa::cuda {

template <std::size_t N>
struct ChandrasekharDevice : FunctionDevice<ChandrasekharDevice<N>, N> {

    double A[N][N];

    __host__ __device__
    void eval_impl(const double* x, double* F) const {
        for (std::size_t i = 0; i < N; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j) s += A[i][j] * x[j];
            F[i] = 1.0 + x[i] * s - x[i];
        }
    }

    __host__ __device__
    void jacobian_impl(const double* x, double* J) const {
        double Ax[N];
        for (std::size_t i = 0; i < N; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j) s += A[i][j] * x[j];
            Ax[i] = s;
        }
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                double v = x[i] * A[i][j];
                if (i == j) v += Ax[i] - 1.0;
                J[i * N + j] = v;
            }
    }

    __host__ __device__
    void hessian_impl(const double*, double* H) const {
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                for (std::size_t k = 0; k < N; ++k) {
                    double v = 0.0;
                    if (i == j) v += A[i][k];
                    if (i == k) v += A[i][j];
                    H[i * N * N + j * N + k] = v;
                }
    }

    // Helper host: inicializar la matriz A de Gauss-Legendre
    void init_A(const double* t, const double* w) {
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                A[i][j] = w[j] * t[i] / (static_cast<double>(N) * (t[i] + t[j]));
    }
};

} // namespace numa::cuda
