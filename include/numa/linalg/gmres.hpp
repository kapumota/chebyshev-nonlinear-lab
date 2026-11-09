#pragma once
// GMRES reiniciado sin matriz con precondicionamiento por la derecha.
// Resuelve A x = b hasta ||b - A x|| <= rtol * ||b||. El conteo de
// productos matriz-vector incluye los usados para residuales de reinicio.

#include "../operators/operator.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace numa::op {

struct GmresResult {
    std::size_t iters = 0;
    std::size_t matvecs = 0;
    std::size_t precs = 0;
    double rel_res = 0.0;
    bool converged = false;
    bool breakdown = false;
};

// Base de Arnoldi del último ciclo de GMRES. La relación guardada es
// A M^{-1} V_k = V_{k+1} Hbar. Sólo es válida mientras A y M no cambien.
struct KrylovSpace {
    std::vector<DVec> V;
    std::vector<std::vector<double>> Hbar;
    std::size_t k = 0;

    void clear() {
        V.clear();
        Hbar.clear();
        k = 0;
    }
};

struct IdentityPrecond {
    void operator()(const DVec& r, DVec& z) const { z = r; }
};

namespace detail {

inline bool valid_krylov_space(const KrylovSpace& space) {
    if (space.k == 0 || space.V.size() < space.k + 1 ||
        space.Hbar.size() < space.k + 1) {
        return false;
    }
    for (std::size_t i = 0; i <= space.k; ++i) {
        if (space.Hbar[i].size() < space.k) return false;
    }
    return true;
}

} // namespace detail

// apply(v, out) calcula A v. precond(r, z) calcula M^{-1} r.
// Si zero_guess es verdadero, la aproximación inicial se fuerza a cero y se
// evita un producto matriz-vector inicial.
// Si save no es nulo, se guarda la base del último ciclo para reciclarla.
// Si r_init no es nulo, se usa como residual inicial y se evita un producto.
template <class Apply, class Precond>
GmresResult gmres(const Apply& apply, const Precond& precond, const DVec& b,
                  DVec& x, double rtol, std::size_t restart,
                  std::size_t max_iters, bool zero_guess = true,
                  KrylovSpace* save = nullptr, const DVec* r_init = nullptr)
{
    const std::size_t n = b.size();
    GmresResult res;
    if (save) save->clear();

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
    if (r_init != nullptr && r_init->size() != n) {
        res.rel_res = 1.0;
        res.breakdown = true;
        return res;
    }
    if (zero_guess) x.assign(n, 0.0);

    DVec r(n), w(n), z(n), tmp(n);
    std::vector<DVec> V(restart + 1, DVec(n));
    std::vector<std::vector<double>> H(
        restart + 1, std::vector<double>(restart, 0.0));
    std::vector<std::vector<double>> Hraw(
        restart + 1, std::vector<double>(restart, 0.0));
    std::vector<double> cs(restart), sn(restart), g(restart + 1), y(restart);

    bool use_given_residual = r_init != nullptr;
    bool first = zero_guess;
    std::size_t k_last = 0;

    while (res.iters < max_iters) {
        if (use_given_residual) {
            r = *r_init;
            use_given_residual = false;
            first = false;
        } else if (first) {
            r = b;
            first = false;
        } else {
            apply(x, w);
            ++res.matvecs;
            for (std::size_t i = 0; i < n; ++i) r[i] = b[i] - w[i];
        }

        const double beta = norm2(r);
        res.rel_res = beta / bnorm;
        if (res.rel_res <= rtol) {
            res.converged = true;
            k_last = 0;
            break;
        }

        for (auto& row : Hraw) std::fill(row.begin(), row.end(), 0.0);
        for (std::size_t i = 0; i < n; ++i) V[0][i] = r[i] / beta;
        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;

        std::size_t k = 0;
        k_last = 0;
        while (k < restart && res.iters < max_iters) {
            const std::size_t j = k;
            precond(V[j], z);
            ++res.precs;
            apply(z, w);
            ++res.matvecs;

            for (std::size_t i = 0; i <= j; ++i) {
                H[i][j] = dot(w, V[i]);
                axpy(-H[i][j], V[i], w);
            }

            H[j + 1][j] = norm2(w);
            if (H[j + 1][j] > 0.0) {
                for (std::size_t i = 0; i < n; ++i) {
                    V[j + 1][i] = w[i] / H[j + 1][j];
                }
            } else {
                std::fill(V[j + 1].begin(), V[j + 1].end(), 0.0);
            }

            for (std::size_t i = 0; i <= j + 1; ++i) {
                Hraw[i][j] = H[i][j];
            }

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
            k_last = k;
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
            for (std::size_t jj = ii + 1; jj < k; ++jj) {
                s -= H[ii][jj] * y[jj];
            }
            y[ii] = s / pivot;
        }

        std::fill(tmp.begin(), tmp.end(), 0.0);
        for (std::size_t i = 0; i < k; ++i) axpy(y[i], V[i], tmp);
        precond(tmp, z);
        ++res.precs;
        axpy(1.0, z, x);

        if (res.rel_res <= rtol) {
            res.converged = true;
            break;
        }
    }

    if (save && k_last > 0) {
        save->k = k_last;
        save->V = V;
        save->Hbar.assign(
            k_last + 1, std::vector<double>(k_last, 0.0));
        for (std::size_t i = 0; i <= k_last; ++i) {
            for (std::size_t j = 0; j < k_last; ++j) {
                save->Hbar[i][j] = Hraw[i][j];
            }
        }
    }

    return res;
}

// Segunda resolución con la misma matriz y otro lado derecho. Proyecta el
// nuevo lado derecho sobre el subespacio guardado sin productos A v y usa
// esa proyección como aproximación inicial. Si el espacio guardado no es
// válido o la proyección es singular, se vuelve a GMRES desde cero.
template <class Apply, class Precond>
GmresResult gmres_recycled(const Apply& apply, const Precond& precond,
                           const DVec& b, DVec& x, double rtol,
                           std::size_t restart, std::size_t max_iters,
                           const KrylovSpace& space,
                           KrylovSpace* save = nullptr)
{
    const std::size_t k = space.k;
    const std::size_t n = b.size();
    x.assign(n, 0.0);

    if (!detail::valid_krylov_space(space)) {
        return gmres(
            apply, precond, b, x, rtol, restart, max_iters, true, save);
    }

    for (std::size_t i = 0; i <= k; ++i) {
        if (space.V[i].size() != n) {
            return gmres(
                apply, precond, b, x, rtol, restart, max_iters, true, save);
        }
    }

    std::vector<double> rhs(k + 1, 0.0);
    for (std::size_t i = 0; i <= k; ++i) {
        rhs[i] = dot(space.V[i], b);
    }

    auto R = space.Hbar;
    std::vector<double> cs(k, 0.0), sn(k, 0.0);
    for (std::size_t j = 0; j < k; ++j) {
        for (std::size_t i = 0; i < j; ++i) {
            const double t = cs[i] * R[i][j] + sn[i] * R[i + 1][j];
            R[i + 1][j] =
                -sn[i] * R[i][j] + cs[i] * R[i + 1][j];
            R[i][j] = t;
        }

        const double denom = std::hypot(R[j][j], R[j + 1][j]);
        if (!(denom > 0.0) || !std::isfinite(denom)) {
            return gmres(
                apply, precond, b, x, rtol, restart, max_iters, true, save);
        }

        cs[j] = R[j][j] / denom;
        sn[j] = R[j + 1][j] / denom;
        R[j][j] = denom;
        R[j + 1][j] = 0.0;

        const double t = cs[j] * rhs[j] + sn[j] * rhs[j + 1];
        rhs[j + 1] = -sn[j] * rhs[j] + cs[j] * rhs[j + 1];
        rhs[j] = t;
    }

    std::vector<double> y(k, 0.0);
    for (std::size_t ii = k; ii-- > 0;) {
        const double pivot = R[ii][ii];
        if (!(std::abs(pivot) > 0.0) || !std::isfinite(pivot)) {
            x.assign(n, 0.0);
            return gmres(
                apply, precond, b, x, rtol, restart, max_iters, true, save);
        }

        double s = rhs[ii];
        for (std::size_t jj = ii + 1; jj < k; ++jj) {
            s -= R[ii][jj] * y[jj];
        }
        y[ii] = s / pivot;
    }

    // r = b - V_{k+1} Hbar y. La proyección no realiza productos A v.
    DVec residual = b;
    DVec preconditioned_basis_combination(n, 0.0);
    DVec projected_x(n, 0.0);

    for (std::size_t i = 0; i <= k; ++i) {
        double coefficient = 0.0;
        for (std::size_t j = 0; j < k; ++j) {
            coefficient += space.Hbar[i][j] * y[j];
        }
        axpy(-coefficient, space.V[i], residual);
    }

    for (std::size_t j = 0; j < k; ++j) {
        axpy(y[j], space.V[j], preconditioned_basis_combination);
    }

    precond(preconditioned_basis_combination, projected_x);
    x = std::move(projected_x);

    auto result = gmres(
        apply, precond, b, x, rtol, restart, max_iters, false, save, &residual);
    ++result.precs;
    return result;
}

} // namespace numa::op
