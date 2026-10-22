#pragma once
// GMRES reiniciado sin matriz con precondicionamiento por la derecha.
// Resuelve A x = b hasta ||b - A x|| <= rtol * ||b||. El conteo de
// productos matriz-vector incluye los usados para residuales de reinicio.

#include "../operators/operator.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace numa::op {

struct GmresResult {
    std::size_t iters = 0;
    std::size_t matvecs = 0;
    double rel_res = 0.0;
    bool converged = false;
    bool breakdown = false;
};

struct IdentityPrecond {
    void operator()(const DVec& r, DVec& z) const { z = r; }
};

// apply(v, out) calcula A v. precond(r, z) calcula M^{-1} r.
// Si zero_guess es verdadero, la aproximación inicial se fuerza a cero y se
// evita un producto matriz-vector inicial.
template <class Apply, class Precond>
GmresResult gmres(const Apply& apply, const Precond& precond, const DVec& b,
                  DVec& x, double rtol, std::size_t restart,
                  std::size_t max_iters, bool zero_guess = true)
{
    const std::size_t n = b.size();
    GmresResult res;
    const double bnorm = norm2(b);
    if (bnorm == 0.0) {
        x.assign(n, 0.0);
        res.converged = true;
        return res;
    }
    if (restart == 0) {
        res.rel_res = 1.0;
        res.breakdown = true;
        return res;
    }
    if (max_iters == 0) {
        res.rel_res = 1.0;
        return res;
    }
    if (zero_guess) x.assign(n, 0.0);

    DVec r(n), w(n), z(n), tmp(n);
    std::vector<DVec> V(restart + 1, DVec(n));
    std::vector<std::vector<double>> H(restart + 1, std::vector<double>(restart, 0.0));
    std::vector<double> cs(restart), sn(restart), g(restart + 1), y(restart);

    bool first = zero_guess;
    while (res.iters < max_iters) {
        if (first) {
            r = b;
            first = false;
        } else {
            apply(x, w);
            ++res.matvecs;
            for (std::size_t i = 0; i < n; ++i) r[i] = b[i] - w[i];
        }
        const double beta = norm2(r);
        res.rel_res = beta / bnorm;
        if (res.rel_res <= rtol) { res.converged = true; break; }

        for (std::size_t i = 0; i < n; ++i) V[0][i] = r[i] / beta;
        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;

        std::size_t k = 0;
        while (k < restart && res.iters < max_iters) {
            const std::size_t j = k;
            precond(V[j], z);
            apply(z, w);
            ++res.matvecs;
            for (std::size_t i = 0; i <= j; ++i) {
                H[i][j] = dot(w, V[i]);
                axpy(-H[i][j], V[i], w);
            }
            H[j + 1][j] = norm2(w);
            if (H[j + 1][j] > 0.0)
                for (std::size_t i = 0; i < n; ++i) V[j + 1][i] = w[i] / H[j + 1][j];

            for (std::size_t i = 0; i < j; ++i) {
                const double t = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                H[i][j] = t;
            }
            const double denom = std::hypot(H[j][j], H[j + 1][j]);
            if (!(denom > 0.0) || !std::isfinite(denom)) {
                res.breakdown = true;
                return res;
            }
            cs[j] = H[j][j] / denom;
            sn[j] = H[j + 1][j] / denom;
            H[j][j] = denom;
            H[j + 1][j] = 0.0;
            g[j + 1] = -sn[j] * g[j];
            g[j] = cs[j] * g[j];

            ++k;
            ++res.iters;
            res.rel_res = std::abs(g[j + 1]) / bnorm;
            if (res.rel_res <= rtol) break;
        }

        for (std::size_t ii = k; ii-- > 0;) {
            const double pivot = H[ii][ii];
            if (!(std::abs(pivot) > 0.0) || !std::isfinite(pivot)) {
                res.breakdown = true;
                return res;
            }
            double s = g[ii];
            for (std::size_t jj = ii + 1; jj < k; ++jj) s -= H[ii][jj] * y[jj];
            y[ii] = s / pivot;
        }
        std::fill(tmp.begin(), tmp.end(), 0.0);
        for (std::size_t i = 0; i < k; ++i) axpy(y[i], V[i], tmp);
        precond(tmp, z);
        axpy(1.0, z, x);

        if (res.rel_res <= rtol) { res.converged = true; break; }
    }
    return res;
}

} // namespace numa::op
