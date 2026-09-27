#pragma once
#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

namespace numa::detail {

// Nodos y pesos de Gauss-Legendre en [0,1], para N arbitrario.
// Algoritmo de Newton-Raphson sobre el polinomio de Legendre P_N,
// equivalente al "gauleg" clasico. Los pesos NO se dividen entre 2,
// para ser consistentes con la tabla de la tesis (N=8), que usa los
// pesos crudos de [-1,1] y solo transforma los nodos a [0,1].
//
// Verificado contra la tabla de la tesis para N=8: coincide en las
// 10 primeras cifras decimales con los valores dados en el Anexo.
template <std::size_t N>
std::pair<std::array<double, N>, std::array<double, N>> gauss_legendre_01()
{
    static_assert(N >= 1, "N debe ser >= 1");
    std::array<double, N> x{}, w{};
    constexpr double pi = 3.14159265358979323846;
    constexpr double eps = 1.0e-15;
    constexpr int max_newton_iter = 100;

    const std::size_t m = (N + 1) / 2; // simetria: solo hay que resolver la mitad

    for (std::size_t i = 0; i < m; ++i) {
        double z = std::cos(pi * (static_cast<double>(i) + 0.75) /
                             (static_cast<double>(N) + 0.5));
        double z1 = 0.0, pp = 0.0;
        int iter = 0;
        do {
            double p1 = 1.0, p2 = 0.0;
            for (std::size_t j = 0; j < N; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * static_cast<double>(j) + 1.0) * z * p2 -
                      static_cast<double>(j) * p3) /
                     (static_cast<double>(j) + 1.0);
            }
            pp = static_cast<double>(N) * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (std::abs(z - z1) > eps && ++iter < max_newton_iter);

        x[i] = -z;
        x[N - 1 - i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
        w[N - 1 - i] = w[i];
    }

    // Transformar nodos de [-1,1] a [0,1]. Los pesos quedan sin dividir
    // entre 2, igual que en la tabla original de la tesis.
    std::array<double, N> t{};
    for (std::size_t i = 0; i < N; ++i) t[i] = (x[i] + 1.0) / 2.0;

    return {t, w};
}

} // namespace numa::detail

namespace numa {

// Uso: GaussLegendre<N>::t().data(), GaussLegendre<N>::w().data()
//
// Nota de migracion respecto a la version original (solo N=8):
// antes `t` y `w` eran arrays estaticos constexpr. Aqui son funciones,
// para poder calcularlos en tiempo de ejecucion para cualquier N sin
// caer en problemas de orden de inicializacion estatica entre
// traducciones. El calculo se hace una sola vez por N gracias al
// "magic static" de C++11+, es seguro entre hilos.
//
// Sitios que deben actualizarse al integrar este archivo:
//   include/numa/problems/chandrasekhar.hpp   (ya actualizado en este paquete)
//   apps/demo_chebyshev_gpu.cu                (ya actualizado en este paquete)
template <std::size_t N>
struct GaussLegendre {
    static const std::array<double, N>& t() {
        static const auto data = detail::gauss_legendre_01<N>();
        return data.first;
    }
    static const std::array<double, N>& w() {
        static const auto data = detail::gauss_legendre_01<N>();
        return data.second;
    }
};

} // namespace numa
