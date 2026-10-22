#pragma once
// Métodos no lineales sin matriz basados en GMRES (R1).
//   newton_krylov       Newton inexacto con términos forzantes de Eisenstat-Walker
//   chebyshev_inexact   Chebyshev con dos resoluciones lineales inexactas
//
// Convención. J z = F(x) y el paso es x - z. En Chebyshev,
//   x+ = x - z - (1/2) u,  con J u = F''(x)[z, z].
// Esta primera entrega no recicla subespacios entre las dos resoluciones.
// Los hiperparámetros son valores iniciales y deben fijarse en el preregistro
// de cada campaña, no ajustarse sobre los resultados.

#include "../linalg/gmres.hpp"
#include "../operators/operator.hpp"
#include <algorithm>
#include <cmath>

namespace numa::op {

struct NoPrecond {
    IdentityPrecond make(const DVec&) const { return {}; }
};

struct KrylovParams {
    double rtol_nl = 1e-8;          // ||F|| <= rtol_nl * ||F(x0)||
    std::size_t max_nl = 50;
    double eta_init = 0.5;          // término forzante inicial
    double eta_max = 0.9;
    double gamma = 0.9;             // Eisenstat-Walker, elección 2
    double alpha = 2.0;
    double eta2_factor = 10.0;      // tolerancia de la segunda resolución
    double eta2_max = 0.5;
    std::size_t gmres_restart = 50;
    std::size_t gmres_max_iters = 500;
    double armijo = 1e-4;
    double min_lambda = 1e-4;
};

struct KrylovResult {
    DVec x;
    bool converged = false;
    std::size_t iters = 0;
    std::size_t n_F = 0, n_JVP = 0, n_D2 = 0;
    std::size_t n_lin_iters = 0;
    std::size_t n_lin_failures = 0;
    double res0 = 0.0, res = 0.0;
};

namespace detail {

// Elección 2 de Eisenstat-Walker con las salvaguardas de Kelley.
inline double forcing(const KrylovParams& p, std::size_t k, double rnorm,
                      double rprev, double eta_prev, double res0)
{
    if (k == 0) return p.eta_init;
    double eta = p.gamma * std::pow(rnorm / rprev, p.alpha);
    const double safe = p.gamma * std::pow(eta_prev, p.alpha);
    if (safe > 0.1) eta = std::max(eta, safe);
    eta = std::min(eta, p.eta_max);
    // Evita sobreresolver más allá de la tolerancia no lineal.
    eta = std::max(eta, 0.5 * p.rtol_nl * res0 / rnorm);
    return std::min(eta, p.eta_max);
}

} // namespace detail

template <NonlinearOperator Op, class PrecondProvider = NoPrecond>
KrylovResult newton_krylov(const Op& f, DVec x, const KrylovParams& p = {},
                           const PrecondProvider& pc = {})
{
    const std::size_t n = f.size();
    KrylovResult out;
    DVec Fx(n), Fn(n), z(n), xn(n);

    f.eval(x, Fx); ++out.n_F;
    double rnorm = norm2(Fx);
    out.res0 = rnorm;
    double rprev = rnorm, eta_prev = p.eta_init;

    for (std::size_t k = 0; k < p.max_nl; ++k) {
        if (rnorm <= p.rtol_nl * out.res0) { out.converged = true; break; }
        const double eta = detail::forcing(p, k, rnorm, rprev, eta_prev, out.res0);
        const auto M = pc.make(x);
        auto apply = [&](const DVec& v, DVec& o) { f.jvp(x, v, o); };

        z.assign(n, 0.0);
        const auto g = gmres(apply, M, Fx, z, eta, p.gmres_restart, p.gmres_max_iters);
        out.n_JVP += g.matvecs;
        out.n_lin_iters += g.iters;
        if (!g.converged) ++out.n_lin_failures;
        if (g.breakdown) break;

        double lambda = 1.0;
        bool ok = false;
        while (lambda >= p.min_lambda) {
            for (std::size_t i = 0; i < n; ++i) xn[i] = x[i] - lambda * z[i];
            f.eval(xn, Fn); ++out.n_F;
            if (norm2(Fn) <= (1.0 - p.armijo * lambda) * rnorm) { ok = true; break; }
            lambda *= 0.5;
        }
        if (!ok) break;

        x = xn; Fx = Fn;
        rprev = rnorm; rnorm = norm2(Fx); eta_prev = eta;
        ++out.iters;
    }
    if (rnorm <= p.rtol_nl * out.res0) out.converged = true;
    out.x = std::move(x);
    out.res = rnorm;
    return out;
}

template <SecondOrderOperator Op, class PrecondProvider = NoPrecond>
KrylovResult chebyshev_inexact(const Op& f, DVec x, const KrylovParams& p = {},
                               const PrecondProvider& pc = {})
{
    const std::size_t n = f.size();
    KrylovResult out;
    DVec Fx(n), Fn(n), z(n), u(n), w(n), xn(n);

    f.eval(x, Fx); ++out.n_F;
    double rnorm = norm2(Fx);
    out.res0 = rnorm;
    double rprev = rnorm, eta_prev = p.eta_init;

    for (std::size_t k = 0; k < p.max_nl; ++k) {
        if (rnorm <= p.rtol_nl * out.res0) { out.converged = true; break; }
        const double eta1 = detail::forcing(p, k, rnorm, rprev, eta_prev, out.res0);
        const double eta2 = std::min(p.eta2_max, p.eta2_factor * eta1);
        const auto M = pc.make(x);
        auto apply = [&](const DVec& v, DVec& o) { f.jvp(x, v, o); };

        z.assign(n, 0.0);
        const auto g1 = gmres(apply, M, Fx, z, eta1, p.gmres_restart, p.gmres_max_iters);
        out.n_JVP += g1.matvecs;
        out.n_lin_iters += g1.iters;
        if (!g1.converged) ++out.n_lin_failures;
        if (g1.breakdown) break;

        f.d2(x, z, w); ++out.n_D2;
        u.assign(n, 0.0);
        const auto g2 = gmres(apply, M, w, u, eta2, p.gmres_restart, p.gmres_max_iters);
        out.n_JVP += g2.matvecs;
        out.n_lin_iters += g2.iters;
        if (!g2.converged) ++out.n_lin_failures;

        // Si la segunda resolución no sufre una ruptura numérica, se intenta primero el
        // paso completo de Chebyshev. En caso contrario, o si ese paso no
        // reduce el residual, se usa el paso de Newton con retroceso.
        bool ok = false;
        if (!g2.breakdown) {
            for (std::size_t i = 0; i < n; ++i) xn[i] = x[i] - z[i] - 0.5 * u[i];
            f.eval(xn, Fn); ++out.n_F;
            if (norm2(Fn) < rnorm) ok = true;
        }
        if (!ok) {
            double lambda = 1.0;
            while (lambda >= p.min_lambda) {
                for (std::size_t i = 0; i < n; ++i) xn[i] = x[i] - lambda * z[i];
                f.eval(xn, Fn); ++out.n_F;
                if (norm2(Fn) <= (1.0 - p.armijo * lambda) * rnorm) { ok = true; break; }
                lambda *= 0.5;
            }
        }
        if (!ok) break;

        x = xn; Fx = Fn;
        rprev = rnorm; rnorm = norm2(Fx); eta_prev = eta1;
        ++out.iters;
    }
    if (rnorm <= p.rtol_nl * out.res0) out.converged = true;
    out.x = std::move(x);
    out.res = rnorm;
    return out;
}

} // namespace numa::op
