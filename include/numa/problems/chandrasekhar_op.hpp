#pragma once

// Ecuación H de Chandrasekhar con dimensión dinámica y cuadratura
// Gauss-Legendre sobre [0,1].
//
// La parametrización usada es:
//
//   x(s) = 1 + (c/2) x(s) integral_0^1 [s/(s+t)] x(t) dt,
//
// con c en [0,1]. El caso histórico de la tesis usa lambda = 1/4,
// equivalente aquí a c = 1/2.
//
// Después de discretizar con Gauss-Legendre:
//
//   F_i(x) = 1 + x_i (A x)_i - x_i,
//   A_ij   = (c/2) w_j mu_i / (mu_i + mu_j).
//
// A no se materializa. Su acción usa memoria O(N) y trabajo O(N^2).
// B5 se emplea principalmente como benchmark controlado de
// condicionamiento, no como evidencia principal de escalabilidad.

#include "../functions/gauss_legendre_dynamic.hpp"
#include "../operators/operator.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace numa::op {

class ChandrasekharOp {
public:
    ChandrasekharOp(std::size_t n, double c)
        : n_(n), c_(c)
    {
        if (n_ == 0) {
            throw std::invalid_argument(
                "ChandrasekharOp requiere N > 0");
        }
        if (!std::isfinite(c_) || c_ < 0.0 || c_ > 1.0) {
            throw std::invalid_argument(
                "ChandrasekharOp requiere c finito en [0,1]");
        }

        const auto rule = quadrature::gauss_legendre_unit(n_);
        mu_ = rule.nodes;
        weights_ = rule.weights;
    }

    std::size_t size() const { return n_; }
    double c() const { return c_; }

    const std::vector<double>& nodes() const { return mu_; }
    const std::vector<double>& weights() const { return weights_; }

    void eval(const DVec& x, DVec& out) const {
        validate(x);
        apply_A(x, Ax_);
        last_x_ = x;
        cache_valid_ = true;

        out.resize(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            out[i] = 1.0 + x[i] * Ax_[i] - x[i];
        }
    }

    void jvp(const DVec& x, const DVec& v, DVec& out) const {
        validate(x);
        validate(v);
        ensure_Ax(x);

        DVec Av;
        apply_A(v, Av);
        out.resize(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            out[i] = v[i] * Ax_[i] + x[i] * Av[i] - v[i];
        }
    }

    void d2(const DVec& x, const DVec& v, DVec& out) const {
        validate(x);
        validate(v);

        DVec Av;
        apply_A(v, Av);
        out.resize(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            out[i] = 2.0 * v[i] * Av[i];
        }
    }

private:
    void validate(const DVec& v) const {
        if (v.size() != n_) {
            throw std::invalid_argument(
                "ChandrasekharOp: dimensión incompatible");
        }
    }

    void ensure_Ax(const DVec& x) const {
        if (!cache_valid_ || last_x_ != x) {
            apply_A(x, Ax_);
            last_x_ = x;
            cache_valid_ = true;
        }
    }

    void apply_A(const DVec& v, DVec& out) const {
        out.assign(n_, 0.0);
        const long n = static_cast<long>(n_);

#pragma omp parallel for schedule(static)
        for (long ii = 0; ii < n; ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            const double mu_i = mu_[i];
            double sum = 0.0;
            for (std::size_t j = 0; j < n_; ++j) {
                sum +=
                    0.5 * c_ * weights_[j] *
                    mu_i / (mu_i + mu_[j]) * v[j];
            }
            out[i] = sum;
        }
    }

    std::size_t n_;
    double c_;
    std::vector<double> mu_;
    std::vector<double> weights_;
    mutable bool cache_valid_ = false;
    mutable DVec last_x_;
    mutable DVec Ax_;
};

} // namespace numa::op
