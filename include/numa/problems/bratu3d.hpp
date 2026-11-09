#pragma once
// Bratu 3D, ecuación de Liouville-Bratu en el cubo unitario con condición
// de Dirichlet nula, discretizada con stencil de 7 puntos.
//   F(u) = -Laplaciano_h(u) - lambda * exp(u)
// N = n^3 incógnitas. Familia B1 del protocolo de gran escala.

#include "../operators/operator.hpp"
#include <cmath>

namespace numa::op {

class Bratu3D {
public:
    explicit Bratu3D(std::size_t n, double lambda = 1.0)
        : n_(n), lambda_(lambda), h2inv_((n + 1.0) * (n + 1.0)) {}

    std::size_t size() const { return n_ * n_ * n_; }
    std::size_t grid() const { return n_; }
    double lambda() const { return lambda_; }

    void eval(const DVec& u, DVec& out) const {
        lap(u, out);
        const long count = static_cast<long>(out.size());
#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (count > 4096)
#endif
        for (long ii = 0; ii < count; ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            out[i] -= lambda_ * std::exp(u[i]);
        }
    }

    void jvp(const DVec& u, const DVec& v, DVec& out) const {
        lap(v, out);
        const long count = static_cast<long>(out.size());
#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (count > 4096)
#endif
        for (long ii = 0; ii < count; ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            out[i] -= lambda_ * std::exp(u[i]) * v[i];
        }
    }

    // F''(u)[v, v] es diagonal: -lambda exp(u) v^2.
    void d2(const DVec& u, const DVec& v, DVec& out) const {
        out.resize(u.size());
        const long count = static_cast<long>(out.size());
#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (count > 4096)
#endif
        for (long ii = 0; ii < count; ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            out[i] = -lambda_ * std::exp(u[i]) * v[i] * v[i];
        }
    }

    // Diagonal del Jacobiano, usada por el precondicionador de Jacobi.
    void jacobi_diag(const DVec& u, DVec& d) const {
        d.resize(u.size());
        const long count = static_cast<long>(d.size());
#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (count > 4096)
#endif
        for (long ii = 0; ii < count; ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            d[i] = 6.0 * h2inv_ - lambda_ * std::exp(u[i]);
        }
    }

private:
    // out = (-Laplaciano_h) v
    void lap(const DVec& v, DVec& out) const {
        out.resize(v.size());
        const std::size_t n = n_;
        const long nl = static_cast<long>(n);
#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for collapse(2) schedule(static) if (nl > 16)
#endif
        for (long ii = 0; ii < nl; ++ii)
            for (long jj = 0; jj < nl; ++jj)
                for (std::size_t k = 0; k < n; ++k) {
                    const std::size_t i = static_cast<std::size_t>(ii);
                    const std::size_t j = static_cast<std::size_t>(jj);
                    const std::size_t id = (i * n + j) * n + k;
                    double s = 6.0 * v[id];
                    if (i > 0)     s -= v[id - n * n];
                    if (i + 1 < n) s -= v[id + n * n];
                    if (j > 0)     s -= v[id - n];
                    if (j + 1 < n) s -= v[id + n];
                    if (k > 0)     s -= v[id - 1];
                    if (k + 1 < n) s -= v[id + 1];
                    out[id] = h2inv_ * s;
                }
    }

    std::size_t n_;
    double lambda_;
    double h2inv_;
};

// Proveedor de precondicionador de Jacobi para Bratu3D. make(x) fija la
// diagonal en el punto x una vez por iteración no lineal.
struct BratuJacobi {
    const Bratu3D& f;

    struct Bound {
        DVec d;
        void operator()(const DVec& r, DVec& z) const {
            z.resize(r.size());
            for (std::size_t i = 0; i < r.size(); ++i) z[i] = r[i] / d[i];
        }
    };

    Bound make(const DVec& x) const {
        Bound b;
        f.jacobi_diag(x, b.d);
        return b;
    }
};

} // namespace numa::op
