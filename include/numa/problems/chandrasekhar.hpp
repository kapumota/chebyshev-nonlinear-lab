#pragma once
#include "../functions/function.hpp"
#include "../functions/gauss_legendre.hpp"

namespace numa::problems {

template <std::size_t N>
struct Chandrasekhar final : IFunction<N> {
    Matrix<N, N> A{};

    Chandrasekhar() {
        // API actualizada: GaussLegendre<N>::t/w ahora son funciones
        // (ver include/numa/functions/gauss_legendre.hpp), porque el
        // calculo para N arbitrario ya no es constexpr.
        const auto& t = GaussLegendre<N>::t();
        const auto& w = GaussLegendre<N>::w();
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                A[i][j] = w[j] * t[i] / (static_cast<double>(N) * (t[i] + t[j]));
    }

    Vec<N> eval(const Vec<N>& x) const override {
        Vec<N> F{};
        for (std::size_t i = 0; i < N; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j) s += A[i][j] * x[j];
            F[i] = 1.0 + x[i] * s - x[i];
        }
        return F;
    }

    Matrix<N, N> jacobian(const Vec<N>& x) const override {
        Vec<N> Ax{};
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                Ax[i] += A[i][j] * x[j];
        Matrix<N, N> J{};
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j)
                J[i][j] = x[i] * A[i][j];
            J[i][i] += Ax[i] - 1.0;
        }
        return J;
    }

    Matrix<N, N * N> hessian(const Vec<N>&) const override {
        Matrix<N, N * N> H{};
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                for (std::size_t k = 0; k < N; ++k) {
                    double v = 0.0;
                    if (i == j) v += A[i][k];
                    if (i == k) v += A[i][j];
                    H[i][j * N + k] = v;
                }
        return H;
    }

    double lipschitz_k() const override { return 0.5; }

    // Acceso por fila, para MPI (Proyecto 3) y para cualquier caso
    // donde N sea grande y calcular F/J/H completos sea el cuello de
    // botella. Costo O(N) por fila, en vez de O(N^2) o O(N^3) de
    // evaluar la funcion completa.
    double row(std::size_t i, const Vec<N>& x) const override {
        double s = 0.0;
        for (std::size_t j = 0; j < N; ++j) s += A[i][j] * x[j];
        return 1.0 + x[i] * s - x[i];
    }

    std::array<double, N> jacobian_row(std::size_t i, const Vec<N>& x) const override {
        double Ax_i = 0.0;
        for (std::size_t j = 0; j < N; ++j) Ax_i += A[i][j] * x[j];
        std::array<double, N> Ji{};
        for (std::size_t j = 0; j < N; ++j) Ji[j] = x[i] * A[i][j];
        Ji[i] += Ax_i - 1.0;
        return Ji;
    }

    std::array<double, N * N> hessian_row(std::size_t i, const Vec<N>&) const override {
        std::array<double, N * N> Hi{};
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t k = 0; k < N; ++k) {
                double v = 0.0;
                if (i == j) v += A[i][k];
                if (i == k) v += A[i][j];
                Hi[j * N + k] = v;
            }
        return Hi;
    }
};

} // namespace numa::problems
