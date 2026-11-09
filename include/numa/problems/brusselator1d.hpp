#pragma once

// Paso implícito de Euler para el Brusselator reacción-difusión 1D usado en
// los tutoriales públicos de PETSc. Este caso es Bext: una aplicación
// canónica no manufacturada. No posee una solución exacta construida para el
// solver.
//
//   u_t - alpha u_xx = A + u^2 v - (B+1) u
//   v_t - alpha v_xx = B u - u^2 v
//
// Se fijan A=1, B=3 y alpha=1/50, con u=1 y v=3 en las fronteras. El estado
// anterior usa u(x)=1+sin(2*pi*x), v(x)=3. El sistema no lineal corresponde a
// un único paso de Euler implícito de tamaño dt.

#include "../operators/operator.hpp"

#include <cmath>
#include <stdexcept>

namespace numa::op {

class Brusselator1DStep {
public:
    explicit Brusselator1DStep(std::size_t n, double dt)
        : n_(n), dt_(dt), h_(1.0 / (static_cast<double>(n) + 1.0)),
          h2inv_(1.0 / (h_ * h_)), previous_(2 * n, 0.0)
    {
        if (n_ == 0 || !(dt_ > 0.0) || !std::isfinite(dt_)) {
            throw std::invalid_argument("Brusselator1DStep requiere n y dt válidos");
        }
        const double pi = std::acos(-1.0);
        for (std::size_t i = 0; i < n_; ++i) {
            const double x = (static_cast<double>(i) + 1.0) * h_;
            previous_[i] = 1.0 + std::sin(2.0 * pi * x);
            previous_[n_ + i] = 3.0;
        }
    }

    std::size_t size() const { return 2 * n_; }
    std::size_t grid() const { return n_; }
    double dt() const { return dt_; }
    static constexpr double A() { return 1.0; }
    static constexpr double B() { return 3.0; }
    static constexpr double diffusion() { return 1.0 / 50.0; }
    const DVec& previous() const { return previous_; }

    void eval(const DVec& w, DVec& out) const {
        validate(w);
        out.assign(size(), 0.0);
        const double invdt = 1.0 / dt_;

#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n_ > 4096)
#endif
        for (long ii = 0; ii < static_cast<long>(n_); ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            const double u = w[i];
            const double v = w[n_ + i];
            const double lap_u = lap_with_boundary(w, 0, i, 1.0);
            const double lap_v = lap_with_boundary(w, 1, i, 3.0);

            out[i] = (u - previous_[i]) * invdt
                   - diffusion() * lap_u
                   - A() - u * u * v + (B() + 1.0) * u;
            out[n_ + i] = (v - previous_[n_ + i]) * invdt
                        - diffusion() * lap_v
                        - B() * u + u * u * v;
        }
    }

    void jvp(const DVec& w, const DVec& d, DVec& out) const {
        validate(w);
        validate(d);
        out.assign(size(), 0.0);
        const double invdt = 1.0 / dt_;

#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n_ > 4096)
#endif
        for (long ii = 0; ii < static_cast<long>(n_); ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            const double u = w[i];
            const double v = w[n_ + i];
            const double du = d[i];
            const double dv = d[n_ + i];
            const double lap_du = lap_direction(d, 0, i);
            const double lap_dv = lap_direction(d, 1, i);

            out[i] = invdt * du
                   - diffusion() * lap_du
                   - (2.0 * u * v) * du - u * u * dv
                   + (B() + 1.0) * du;
            out[n_ + i] = invdt * dv
                        - diffusion() * lap_dv
                        - B() * du
                        + (2.0 * u * v) * du + u * u * dv;
        }
    }

    void d2(const DVec& w, const DVec& d, DVec& out) const {
        validate(w);
        validate(d);
        out.assign(size(), 0.0);

#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n_ > 4096)
#endif
        for (long ii = 0; ii < static_cast<long>(n_); ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            const double u = w[i];
            const double v = w[n_ + i];
            const double du = d[i];
            const double dv = d[n_ + i];
            const double q = 2.0 * v * du * du + 4.0 * u * du * dv;
            out[i] = -q;
            out[n_ + i] = q;
        }
    }

    void local_jacobi_blocks(
        const DVec& w,
        DVec& a11,
        DVec& a12,
        DVec& a21,
        DVec& a22) const
    {
        validate(w);
        a11.resize(n_);
        a12.resize(n_);
        a21.resize(n_);
        a22.resize(n_);
        const double base = 1.0 / dt_ + 2.0 * diffusion() * h2inv_;
        for (std::size_t i = 0; i < n_; ++i) {
            const double u = w[i];
            const double v = w[n_ + i];
            a11[i] = base - 2.0 * u * v + (B() + 1.0);
            a12[i] = -u * u;
            a21[i] = -B() + 2.0 * u * v;
            a22[i] = base + u * u;
        }
    }

private:
    void validate(const DVec& w) const {
        if (w.size() != size()) {
            throw std::invalid_argument("Brusselator1DStep: dimensión incompatible");
        }
    }

    double value_at(const DVec& w, std::size_t comp, long i) const {
        if (i < 0 || i >= static_cast<long>(n_)) {
            return comp == 0 ? 1.0 : 3.0;
        }
        return w[comp * n_ + static_cast<std::size_t>(i)];
    }

    double direction_component(const DVec& d, std::size_t component, long i) const {
        if (i < 0 || i >= static_cast<long>(n_)) return 0.0;
        return d[component * n_ + static_cast<std::size_t>(i)];
    }

    double lap_with_boundary(
        const DVec& w,
        std::size_t component,
        std::size_t i,
        double) const
    {
        const long ii = static_cast<long>(i);
        const double center = w[component * n_ + i];
        return (value_at(w, component, ii - 1)
              - 2.0 * center
              + value_at(w, component, ii + 1)) * h2inv_;
    }

    double lap_direction(
        const DVec& d,
        std::size_t component,
        std::size_t i) const
    {
        const long ii = static_cast<long>(i);
        const double center = d[component * n_ + i];
        return (direction_component(d, component, ii - 1)
              - 2.0 * center
              + direction_component(d, component, ii + 1)) * h2inv_;
    }

    std::size_t n_;
    double dt_;
    double h_;
    double h2inv_;
    DVec previous_;
};

struct BrusselatorJacobi {
    const Brusselator1DStep& f;

    struct Bound {
        DVec a11, a12, a21, a22;

        void operator()(const DVec& r, DVec& z) const {
            const std::size_t n = a11.size();
            z.assign(2 * n, 0.0);
            for (std::size_t i = 0; i < n; ++i) {
                const double det = a11[i] * a22[i] - a12[i] * a21[i];
                if (!std::isfinite(det) || std::abs(det) < 1e-14) {
                    z[i] = r[i];
                    z[n + i] = r[n + i];
                    continue;
                }
                z[i] = (a22[i] * r[i] - a12[i] * r[n + i]) / det;
                z[n + i] = (-a21[i] * r[i] + a11[i] * r[n + i]) / det;
            }
        }
    };

    Bound make(const DVec& x) const {
        Bound bound;
        f.local_jacobi_blocks(
            x, bound.a11, bound.a12, bound.a21, bound.a22);
        return bound;
    }
};

} // namespace numa::op
