#pragma once

// Cuadratura de Gauss-Legendre dinámica sobre [0,1].
//
// El baseline histórico conserva GaussLegendre<8> con la tabla de la tesis.
// Esta implementación genera nodos y pesos para N arbitrario sin alterar ese
// baseline. Los ceros de P_N se obtienen con Newton y se explota la simetría.

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace numa::quadrature {

struct GaussLegendreRule {
    std::vector<double> nodes;
    std::vector<double> weights;
};

namespace detail {

inline std::pair<long double, long double>
legendre_value_derivative(std::size_t n, long double x)
{
    long double p0 = 1.0L;
    if (n == 0) return {p0, 0.0L};

    long double p1 = x;
    if (n == 1) return {p1, 1.0L};

    for (std::size_t k = 2; k <= n; ++k) {
        const long double kd = static_cast<long double>(k);
        const long double pk =
            ((2.0L * kd - 1.0L) * x * p1 - (kd - 1.0L) * p0) / kd;
        p0 = p1;
        p1 = pk;
    }

    const long double nd = static_cast<long double>(n);
    const long double derivative =
        nd * (x * p1 - p0) / (x * x - 1.0L);
    return {p1, derivative};
}

} // namespace detail

inline GaussLegendreRule gauss_legendre_unit(std::size_t n)
{
    if (n == 0) {
        throw std::invalid_argument("gauss_legendre_unit requiere N > 0");
    }

    GaussLegendreRule rule;
    rule.nodes.resize(n);
    rule.weights.resize(n);

    const std::size_t half = (n + 1) / 2;
    const long double pi = std::acos(-1.0L);
    const long double eps =
        32.0L * std::numeric_limits<long double>::epsilon();

    for (std::size_t i = 0; i < half; ++i) {
        const long double guess_index =
            static_cast<long double>(i) + 0.75L;
        const long double nd = static_cast<long double>(n);
        long double z =
            std::cos(pi * guess_index / (nd + 0.5L));

        bool converged = false;
        for (std::size_t iter = 0; iter < 100; ++iter) {
            const auto [value, derivative] =
                detail::legendre_value_derivative(n, z);
            if (!std::isfinite(static_cast<double>(value)) ||
                !std::isfinite(static_cast<double>(derivative)) ||
                derivative == 0.0L) {
                throw std::runtime_error(
                    "Fallo numérico al generar Gauss-Legendre");
            }

            const long double delta = value / derivative;
            z -= delta;
            if (std::abs(delta) <= eps * (1.0L + std::abs(z))) {
                converged = true;
                break;
            }
        }

        if (!converged) {
            throw std::runtime_error(
                "Gauss-Legendre no convergió en 100 iteraciones");
        }

        const auto [value, derivative] =
            detail::legendre_value_derivative(n, z);
        (void)value;
        const long double weight_minus_one_one =
            2.0L / ((1.0L - z * z) * derivative * derivative);

        const double left = static_cast<double>(0.5L * (1.0L - z));
        const double right = static_cast<double>(0.5L * (1.0L + z));
        const double weight =
            static_cast<double>(0.5L * weight_minus_one_one);

        rule.nodes[i] = left;
        rule.nodes[n - 1 - i] = right;
        rule.weights[i] = weight;
        rule.weights[n - 1 - i] = weight;
    }

    return rule;
}

} // namespace numa::quadrature
