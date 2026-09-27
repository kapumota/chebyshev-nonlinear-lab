#pragma once
#include "newton.hpp"
#include "../recurrence/sigma.hpp"
#include "../recurrence/pol.hpp"

namespace numa {

// Algoritmo 5: Método iterativo de Chebyshev (secuencial).
// x_{n+1} = x_n - [I + (1/2) Γ_n F''(x_n) Γ_n F(x_n)] Γ_n F(x_n)
// Implementado como:
//   z = Γ_n F(x_n)
//   w[i] = Σ_{j,k} H[i][j*N+k] z[j] z[k]
//   u = Γ_n w
//   dx = z + 0.5 * u
template <std::size_t N>
SolverResult<N> chebyshev(const IFunction<N>& f,
                          Vec<N> x,
                          double tol = 1e-12,
                          std::size_t max_iter = 20)
{
    double er = 0.0;
    std::size_t k = 0;

    while (k < max_iter) {
        const Vec<N> previous_x = x;
        auto F = f.eval(x);
        auto J = f.jacobian(x);
        auto H = f.hessian(x);

        auto Dz = J; auto dz = F;
        gauss_seidel(Dz, dz);
        Vec<N> z = tri_sup(Dz, dz);

        Vec<N> w{};
        for (std::size_t i = 0; i < N; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j)
                for (std::size_t kk = 0; kk < N; ++kk)
                    s += H[i][j * N + kk] * z[j] * z[kk];
            w[i] = s;
        }

        auto Du = J; auto du = w;
        gauss_seidel(Du, du);
        Vec<N> u = tri_sup(Du, du);

        Vec<N> dx{};
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
