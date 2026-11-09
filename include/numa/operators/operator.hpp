#pragma once
// Interfaz de operadores no lineales de dimensión dinámica (R1).
//
// Un operador expone F, el producto Jacobiano-vector y, de forma opcional,
// la segunda derivada direccional F''(x)[v, v]. No materializa Jacobianos
// ni Hessianos. El baseline de tamaño fijo se conserva sin cambios y se
// conecta mediante BaselineAdapter.

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <vector>

namespace numa::op {

using DVec = std::vector<double>;

inline double dot(const DVec& a, const DVec& b) {
    double s = 0.0;
    const long n = static_cast<long>(a.size());
#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for reduction(+ : s) schedule(static) if (n > 4096)
#endif
    for (long i = 0; i < n; ++i) s += a[static_cast<std::size_t>(i)] * b[static_cast<std::size_t>(i)];
    return s;
}

inline double norm2(const DVec& a) { return std::sqrt(dot(a, a)); }

// y += a * x
inline void axpy(double a, const DVec& x, DVec& y) {
    const long n = static_cast<long>(x.size());
#ifdef NUMA_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n > 4096)
#endif
    for (long i = 0; i < n; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        y[index] += a * x[index];
    }
}

// Requisitos mínimos. El llamador dimensiona out con size() elementos.
template <class Op>
concept NonlinearOperator =
    requires(const Op& f, const DVec& x, const DVec& v, DVec& out) {
        { f.size() } -> std::convertible_to<std::size_t>;
        f.eval(x, out);      // out = F(x)
        f.jvp(x, v, out);    // out = F'(x) v
    };

// Operador con segunda derivada direccional, out = F''(x)[v, v].
template <class Op>
concept SecondOrderOperator =
    NonlinearOperator<Op> &&
    requires(const Op& f, const DVec& x, const DVec& v, DVec& out) {
        f.d2(x, v, out);
    };

// Contadores de costo para la contabilidad del protocolo de gran escala.
struct OpCounters {
    std::size_t n_F = 0;
    std::size_t n_JVP = 0;
    std::size_t n_D2 = 0;
};

// Envoltorio que cuenta llamadas sin alterar resultados.
template <class Op>
class Counted {
public:
    explicit Counted(const Op& op) : op_(op) {}

    std::size_t size() const { return op_.size(); }

    void eval(const DVec& x, DVec& out) const {
        ++c_.n_F;
        op_.eval(x, out);
    }

    void jvp(const DVec& x, const DVec& v, DVec& out) const {
        ++c_.n_JVP;
        op_.jvp(x, v, out);
    }

    void d2(const DVec& x, const DVec& v, DVec& out) const
        requires SecondOrderOperator<Op>
    {
        ++c_.n_D2;
        op_.d2(x, v, out);
    }

    const OpCounters& counters() const { return c_; }
    void reset() const { c_ = OpCounters{}; }

private:
    const Op& op_;
    mutable OpCounters c_;
};

// Producto Jacobiano-vector por diferencia hacia adelante. Usa 2 evaluaciones de F.
template <class Op>
void jvp_fd(const Op& f, const DVec& x, const DVec& v, DVec& out) {
    const std::size_t n = x.size();
    out.assign(n, 0.0);
    const double nv = norm2(v);
    if (nv == 0.0) return;
    const double h = std::sqrt(std::numeric_limits<double>::epsilon()) *
                     (1.0 + norm2(x)) / nv;
    DVec xp = x, f0(n), f1(n);
    axpy(h, v, xp);
    f.eval(x, f0);
    f.eval(xp, f1);
    for (std::size_t i = 0; i < n; ++i) out[i] = (f1[i] - f0[i]) / h;
}

// Segunda derivada direccional F''(x)[v, v] por diferencia central. Usa 3
// evaluaciones de F. El paso sigue la escala epsilon^(1/4), a calibrar en C2.
template <class Op>
void d2_fd(const Op& f, const DVec& x, const DVec& v, DVec& out) {
    const std::size_t n = x.size();
    out.assign(n, 0.0);
    const double nv = norm2(v);
    if (nv == 0.0) return;
    const double h = std::pow(std::numeric_limits<double>::epsilon(), 0.25) *
                     (1.0 + norm2(x)) / nv;
    DVec xp = x, xm = x, f0(n), fp(n), fm(n);
    axpy(h, v, xp);
    axpy(-h, v, xm);
    f.eval(x, f0);
    f.eval(xp, fp);
    f.eval(xm, fm);
    for (std::size_t i = 0; i < n; ++i)
        out[i] = (fp[i] - 2.0 * f0[i] + fm[i]) / (h * h);
}

} // namespace numa::op
