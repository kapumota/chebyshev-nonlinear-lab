#pragma once

// Burgers 2D viscoso en un paso de Euler implícito con solución manufacturada.
// B3a es un benchmark controlado para verificar un Jacobiano no simétrico y
// resoluciones lineales más difíciles. No constituye validación externa.
//
//   F(w) = w / dt + C(w) - nu Laplaciano_h(w) - g
//   C(w) = (u u_x + v u_y, u v_x + v v_y)
//
// La fuente g se construye de forma que w* sea una raíz discreta exacta.
// N = 2 n^2. La segunda derivada direccional es exacta porque F es cuadrática.

#include "../linalg/helmholtz_mg2d.hpp"
#include "../operators/operator.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace numa::op {

class Burgers2D {
public:
    Burgers2D(std::size_t n, double nu, double dt, double amplitude = 1.0)
        : n_(n), nu_(nu), dt_(dt), amplitude_(amplitude),
          source_(2 * n * n, 0.0), exact_(2 * n * n, 0.0)
    {
        if (n_ == 0 || !(nu_ > 0.0) || !(dt_ > 0.0) ||
            !std::isfinite(nu_) || !std::isfinite(dt_) ||
            !std::isfinite(amplitude_)) {
            throw std::invalid_argument("Burgers2D requiere n, nu y dt válidos");
        }

        const double pi = std::acos(-1.0);
        const double h = 1.0 / (static_cast<double>(n_) + 1.0);
        const std::size_t nn = n_ * n_;
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t j = 0; j < n_; ++j) {
                const double y = (static_cast<double>(i) + 1.0) * h;
                const double x = (static_cast<double>(j) + 1.0) * h;
                const std::size_t id = i * n_ + j;
                exact_[id] = amplitude_ * std::sin(pi * x) * std::sin(pi * y);
                exact_[nn + id] =
                    amplitude_ * std::sin(2.0 * pi * x) * std::sin(pi * y);
            }
        }

        DVec residual_without_source;
        eval_without_source(exact_, residual_without_source);
        source_ = std::move(residual_without_source);
    }

    std::size_t size() const { return 2 * n_ * n_; }
    std::size_t grid() const { return n_; }
    double nu() const { return nu_; }
    double dt() const { return dt_; }
    double amplitude() const { return amplitude_; }
    double cell_reynolds() const {
        return std::abs(amplitude_) /
               (nu_ * (static_cast<double>(n_) + 1.0));
    }
    const DVec& exact() const { return exact_; }

    void eval(const DVec& w, DVec& out) const {
        validate(w);
        eval_without_source(w, out);
        for (std::size_t i = 0; i < out.size(); ++i) out[i] -= source_[i];
    }

    void jvp(const DVec& w, const DVec& d, DVec& out) const {
        validate(w);
        validate(d);
        const std::size_t nn = n_ * n_;
        const double h = 1.0 / (static_cast<double>(n_) + 1.0);
        const double inv2h = 1.0 / (2.0 * h);
        const double invh2 = 1.0 / (h * h);
        const double invdt = 1.0 / dt_;
        out.assign(size(), 0.0);

        auto at = [&](const DVec& a, std::size_t component, long i, long j) {
            if (i < 0 || j < 0 || i >= static_cast<long>(n_) ||
                j >= static_cast<long>(n_)) {
                return 0.0;
            }
            return a[component * nn + static_cast<std::size_t>(i) * n_ +
                     static_cast<std::size_t>(j)];
        };

#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nn > 4096)
#endif
        for (long ii = 0; ii < static_cast<long>(n_); ++ii) {
            for (long jj = 0; jj < static_cast<long>(n_); ++jj) {
                const std::size_t i = static_cast<std::size_t>(ii);
                const std::size_t j = static_cast<std::size_t>(jj);
                const std::size_t id = i * n_ + j;
                const double u = w[id];
                const double v = w[nn + id];
                const double du = d[id];
                const double dv = d[nn + id];

                const double ux = (at(w, 0, ii, jj + 1) - at(w, 0, ii, jj - 1)) * inv2h;
                const double uy = (at(w, 0, ii + 1, jj) - at(w, 0, ii - 1, jj)) * inv2h;
                const double vx = (at(w, 1, ii, jj + 1) - at(w, 1, ii, jj - 1)) * inv2h;
                const double vy = (at(w, 1, ii + 1, jj) - at(w, 1, ii - 1, jj)) * inv2h;
                const double dux = (at(d, 0, ii, jj + 1) - at(d, 0, ii, jj - 1)) * inv2h;
                const double duy = (at(d, 0, ii + 1, jj) - at(d, 0, ii - 1, jj)) * inv2h;
                const double dvx = (at(d, 1, ii, jj + 1) - at(d, 1, ii, jj - 1)) * inv2h;
                const double dvy = (at(d, 1, ii + 1, jj) - at(d, 1, ii - 1, jj)) * inv2h;

                const double lap_du =
                    (at(d, 0, ii, jj + 1) + at(d, 0, ii, jj - 1) +
                     at(d, 0, ii + 1, jj) + at(d, 0, ii - 1, jj) - 4.0 * du) * invh2;
                const double lap_dv =
                    (at(d, 1, ii, jj + 1) + at(d, 1, ii, jj - 1) +
                     at(d, 1, ii + 1, jj) + at(d, 1, ii - 1, jj) - 4.0 * dv) * invh2;

                out[id] = invdt * du + du * ux + u * dux + dv * uy + v * duy - nu_ * lap_du;
                out[nn + id] = invdt * dv + du * vx + u * dvx + dv * vy + v * dvy - nu_ * lap_dv;
            }
        }
    }

    void d2(const DVec& w, const DVec& d, DVec& out) const {
        validate(w);
        validate(d);
        const std::size_t nn = n_ * n_;
        const double h = 1.0 / (static_cast<double>(n_) + 1.0);
        const double inv2h = 1.0 / (2.0 * h);
        out.assign(size(), 0.0);

        auto at = [&](std::size_t component, long i, long j) {
            if (i < 0 || j < 0 || i >= static_cast<long>(n_) ||
                j >= static_cast<long>(n_)) {
                return 0.0;
            }
            return d[component * nn + static_cast<std::size_t>(i) * n_ +
                     static_cast<std::size_t>(j)];
        };

#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nn > 4096)
#endif
        for (long ii = 0; ii < static_cast<long>(n_); ++ii) {
            for (long jj = 0; jj < static_cast<long>(n_); ++jj) {
                const std::size_t i = static_cast<std::size_t>(ii);
                const std::size_t j = static_cast<std::size_t>(jj);
                const std::size_t id = i * n_ + j;
                const double du = d[id];
                const double dv = d[nn + id];
                const double dux = (at(0, ii, jj + 1) - at(0, ii, jj - 1)) * inv2h;
                const double duy = (at(0, ii + 1, jj) - at(0, ii - 1, jj)) * inv2h;
                const double dvx = (at(1, ii, jj + 1) - at(1, ii, jj - 1)) * inv2h;
                const double dvy = (at(1, ii + 1, jj) - at(1, ii - 1, jj)) * inv2h;
                out[id] = 2.0 * (du * dux + dv * duy);
                out[nn + id] = 2.0 * (du * dvx + dv * dvy);
            }
        }
    }

private:
    void validate(const DVec& w) const {
        if (w.size() != size()) {
            throw std::invalid_argument("Burgers2D: dimensión incompatible");
        }
    }

    void eval_without_source(const DVec& w, DVec& out) const {
        const std::size_t nn = n_ * n_;
        const double h = 1.0 / (static_cast<double>(n_) + 1.0);
        const double inv2h = 1.0 / (2.0 * h);
        const double invh2 = 1.0 / (h * h);
        const double invdt = 1.0 / dt_;
        out.assign(size(), 0.0);

        auto at = [&](std::size_t component, long i, long j) {
            if (i < 0 || j < 0 || i >= static_cast<long>(n_) ||
                j >= static_cast<long>(n_)) {
                return 0.0;
            }
            return w[component * nn + static_cast<std::size_t>(i) * n_ +
                     static_cast<std::size_t>(j)];
        };

#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nn > 4096)
#endif
        for (long ii = 0; ii < static_cast<long>(n_); ++ii) {
            for (long jj = 0; jj < static_cast<long>(n_); ++jj) {
                const std::size_t i = static_cast<std::size_t>(ii);
                const std::size_t j = static_cast<std::size_t>(jj);
                const std::size_t id = i * n_ + j;
                const double u = w[id];
                const double v = w[nn + id];
                const double ux = (at(0, ii, jj + 1) - at(0, ii, jj - 1)) * inv2h;
                const double uy = (at(0, ii + 1, jj) - at(0, ii - 1, jj)) * inv2h;
                const double vx = (at(1, ii, jj + 1) - at(1, ii, jj - 1)) * inv2h;
                const double vy = (at(1, ii + 1, jj) - at(1, ii - 1, jj)) * inv2h;
                const double lap_u =
                    (at(0, ii, jj + 1) + at(0, ii, jj - 1) +
                     at(0, ii + 1, jj) + at(0, ii - 1, jj) - 4.0 * u) * invh2;
                const double lap_v =
                    (at(1, ii, jj + 1) + at(1, ii, jj - 1) +
                     at(1, ii + 1, jj) + at(1, ii - 1, jj) - 4.0 * v) * invh2;
                out[id] = invdt * u + u * ux + v * uy - nu_ * lap_u;
                out[nn + id] = invdt * v + u * vx + v * vy - nu_ * lap_v;
            }
        }
    }

    std::size_t n_;
    double nu_;
    double dt_;
    double amplitude_;
    DVec source_;
    DVec exact_;
};

struct BurgersMGProvider {
    std::shared_ptr<HelmholtzMG2D> mg;
    std::size_t nn = 0;

    BurgersMGProvider(std::size_t n, double nu, double dt)
        : mg(std::make_shared<HelmholtzMG2D>(n, 1.0 / dt, nu)),
          nn(n * n) {}

    struct Bound {
        const HelmholtzMG2D* mg = nullptr;
        std::size_t nn = 0;

        void operator()(const DVec& r, DVec& z) const {
            DVec ru(r.begin(), r.begin() + static_cast<std::ptrdiff_t>(nn));
            DVec rv(r.begin() + static_cast<std::ptrdiff_t>(nn), r.end());
            DVec zu, zv;
            (*mg)(ru, zu);
            (*mg)(rv, zv);
            z.resize(r.size());
            std::copy(zu.begin(), zu.end(), z.begin());
            std::copy(zv.begin(), zv.end(),
                      z.begin() + static_cast<std::ptrdiff_t>(nn));
        }
    };

    Bound make(const DVec&) const { return Bound{mg.get(), nn}; }
};

} // namespace numa::op
