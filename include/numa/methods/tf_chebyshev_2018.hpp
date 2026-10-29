#pragma once

// Baseline estructural del método inexacto tensor-free Chebyshev-Halley de
// Eustaquio, Ribeiro y Dumett (2018), especializado a alpha = 0
// (Chebyshev tensor-free). No se presenta como método nuevo.
//
// Convención local. El primer GMRES resuelve J z = F, por lo que d1 = -z.
// La curvatura se aproxima mediante
//   C(x)d1 = [J(x + h d1)d1 - J(x)d1] / h.
// El segundo sistema, para alpha = 0, es
//   J d2 = -0.5 C(x)d1,
// y la actualización es x+ = x + d1 + d2.

#include "../linalg/gmres.hpp"
#include "krylov_newton.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace numa::op {

inline double tf2018_step_length(const DVec& x, const DVec& d1) {
    const double eps = std::numeric_limits<double>::epsilon();
    const double nd = norm2(d1);
    if (!(nd > 0.0)) return 2.0 * eps;
    const double raw = std::sqrt(eps) * norm2(x) / nd;
    return std::clamp(raw, 2.0 * eps, 0.5);
}

template <NonlinearOperator Op>
void tensor_free_curvature_2018(
    const Op& f,
    const DVec& x,
    const DVec& d1,
    DVec& out,
    double* h_used = nullptr)
{
    const double h = tf2018_step_length(x, d1);
    DVec xh = x;
    axpy(h, d1, xh);
    DVec j0(x.size()), j1(x.size());
    f.jvp(x, d1, j0);
    f.jvp(xh, d1, j1);
    out.resize(x.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = (j1[i] - j0[i]) / h;
    }
    if (h_used != nullptr) *h_used = h;
}

template <NonlinearOperator Op, class PrecondProvider = NoPrecond>
KrylovResult tf_chebyshev_2018(
    const Op& f,
    DVec x,
    const KrylovParams& p = {},
    const PrecondProvider& pc = {})
{
    const std::size_t n = f.size();
    KrylovResult out;
    DVec Fx(n), Fn(n), z(n), d1(n), d2(n), curvature(n), rhs2(n);

    f.eval(x, Fx);
    ++out.n_F;
    double rnorm = norm2(Fx);
    out.res0 = rnorm;
    const double eps = std::numeric_limits<double>::epsilon();

    for (std::size_t k = 0; k < p.max_nl; ++k) {
        if (rnorm <= p.rtol_nl * out.res0) {
            out.converged = true;
            break;
        }

        const double eta1 = std::clamp(
            std::min(1e-8, rnorm * rnorm), 2.0 * eps, 1e-8);
        const double eta2 = std::clamp(
            std::min(1e-8, rnorm), 2.0 * eps, 1e-8);
        const auto M = pc.make(x);
        auto apply = [&](const DVec& v, DVec& o) { f.jvp(x, v, o); };

        z.assign(n, 0.0);
        const auto g1 = gmres(
            apply, M, Fx, z, eta1,
            p.gmres_restart, p.gmres_max_iters);
        out.n_JVP += g1.matvecs;
        out.n_lin_iters += g1.iters;
        out.n_prec += g1.precs;
        if (!g1.converged) {
            ++out.n_lin_failures;
            break;
        }
        if (g1.breakdown) break;

        d1 = z;
        for (double& value : d1) value = -value;

        tensor_free_curvature_2018(f, x, d1, curvature);
        out.n_JVP += 2;
        for (std::size_t i = 0; i < n; ++i) {
            rhs2[i] = -0.5 * curvature[i];
        }

        d2.assign(n, 0.0);
        const auto g2 = gmres(
            apply, M, rhs2, d2, eta2,
            p.gmres_restart, p.gmres_max_iters);
        out.n_JVP += g2.matvecs;
        out.n_lin_iters += g2.iters;
        out.n_prec += g2.precs;
        if (!g2.converged) {
            ++out.n_lin_failures;
            break;
        }
        if (g2.breakdown) break;

        axpy(1.0, d1, x);
        axpy(1.0, d2, x);
        f.eval(x, Fn);
        ++out.n_F;
        Fx = Fn;
        rnorm = norm2(Fx);
        ++out.iters;

        if (!std::isfinite(rnorm)) break;
    }

    if (rnorm <= p.rtol_nl * out.res0) out.converged = true;
    out.x = std::move(x);
    out.res = rnorm;
    return out;
}

} // namespace numa::op
