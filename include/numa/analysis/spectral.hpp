#pragma once

// Utilidades espectrales para C1. La rutina usa iteración inversa con GMRES
// sobre un operador simétrico y devuelve un estimador del autovalor de menor
// magnitud junto con el residual del par propio.

#include "../linalg/gmres.hpp"
#include "../operators/operator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace numa::analysis {

using numa::op::DVec;
using numa::op::dot;
using numa::op::norm2;
using numa::op::axpy;

struct EigenEstimate {
    bool converged = false;
    double eigenvalue = std::numeric_limits<double>::quiet_NaN();
    double residual_rel = std::numeric_limits<double>::infinity();
    std::size_t outer_iterations = 0;
    std::size_t linear_iterations = 0;
    std::size_t linear_matvecs = 0;
    std::size_t preconditioner_applications = 0;
    std::size_t linear_failures = 0;
    DVec vector;
};

template <class Apply, class Precond>
EigenEstimate inverse_smallest_eigenpair(
    const Apply& apply,
    const Precond& precond,
    std::size_t n,
    double eigen_rtol = 1e-7,
    double linear_rtol = 1e-10,
    std::size_t max_outer = 30,
    std::size_t restart = 100,
    std::size_t max_linear = 2000)
{
    EigenEstimate result;
    if (n == 0 || !(eigen_rtol > 0.0) || !(linear_rtol > 0.0) ||
        max_outer == 0 || restart == 0 || max_linear == 0) {
        return result;
    }

    DVec q(n);
    for (std::size_t i = 0; i < n; ++i) {
        q[i] = 1.0 + 0.1 * std::sin(
            0.6180339887498948 * static_cast<double>(i + 1));
    }

    const double initial_norm = norm2(q);
    if (!(initial_norm > 0.0) || !std::isfinite(initial_norm)) {
        return result;
    }
    for (double& value : q) value /= initial_norm;

    double previous = std::numeric_limits<double>::quiet_NaN();
    DVec y(n), jq(n), residual(n);

    for (std::size_t outer = 0; outer < max_outer; ++outer) {
        y.assign(n, 0.0);
        const auto linear = numa::op::gmres(
            apply, precond, q, y, linear_rtol, restart, max_linear);

        result.linear_iterations += linear.iters;
        result.linear_matvecs += linear.matvecs;
        result.preconditioner_applications += linear.precs;
        if (!linear.converged || linear.breakdown) {
            ++result.linear_failures;
            return result;
        }

        const double y_norm = norm2(y);
        if (!(y_norm > std::numeric_limits<double>::min()) ||
            !std::isfinite(y_norm)) {
            ++result.linear_failures;
            return result;
        }
        for (std::size_t i = 0; i < n; ++i) q[i] = y[i] / y_norm;

        apply(q, jq);
        const double theta = dot(q, jq);
        if (!std::isfinite(theta)) return result;

        residual = jq;
        axpy(-theta, q, residual);
        const double denominator = std::max(1.0, std::abs(theta));
        const double residual_rel = norm2(residual) / denominator;

        result.eigenvalue = theta;
        result.residual_rel = residual_rel;
        result.outer_iterations = outer + 1;
        result.vector = q;

        const double change =
            std::isfinite(previous)
                ? std::abs(theta - previous) / denominator
                : std::numeric_limits<double>::infinity();

        if (residual_rel <= eigen_rtol &&
            (!std::isfinite(previous) || change <= 10.0 * eigen_rtol)) {
            result.converged = true;
            return result;
        }
        previous = theta;
    }

    return result;
}

} // namespace numa::analysis
