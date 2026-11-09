#pragma once
// Adaptador de IFunction<N> del baseline a la interfaz de operadores.
// Permite ejecutar los métodos de R1 sobre los casos pequeños y usar el
// baseline como suite de regresión.

#include "operator.hpp"
#include "../functions/function.hpp"

namespace numa::op {

template <std::size_t N>
class BaselineAdapter {
public:
    explicit BaselineAdapter(const IFunction<N>& f) : f_(f) {}

    std::size_t size() const { return N; }

    void eval(const DVec& x, DVec& out) const {
        const auto F = f_.eval(to_vec(x));
        out.assign(F.begin(), F.end());
    }

    void jvp(const DVec& x, const DVec& v, DVec& out) const {
        const auto J = f_.jacobian(to_vec(x));
        out.assign(N, 0.0);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) out[i] += J[i][j] * v[j];
    }

    void d2(const DVec& x, const DVec& v, DVec& out) const {
        const auto H = f_.hessian(to_vec(x));
        out.assign(N, 0.0);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                for (std::size_t k = 0; k < N; ++k)
                    out[i] += H[i][j * N + k] * v[j] * v[k];
    }

private:
    static Vec<N> to_vec(const DVec& x) {
        Vec<N> r{};
        for (std::size_t i = 0; i < N; ++i) r[i] = x[i];
        return r;
    }

    const IFunction<N>& f_;
};

} // namespace numa::op
