#pragma once
#include <cstddef>

// CRTP: polimorfismo en tiempo de compilación para que las funciones
// sean invocables desde kernels __global__ / __device__.
// Las funciones virtuales NO funcionan en device code, por eso CRTP.

namespace numa::cuda {

template <typename Derived, std::size_t N>
struct FunctionDevice {

    __host__ __device__
    void eval(const double* x, double* F) const {
        static_cast<const Derived*>(this)->eval_impl(x, F);
    }

    __host__ __device__
    void jacobian(const double* x, double* J) const {
        static_cast<const Derived*>(this)->jacobian_impl(x, J);
    }

    __host__ __device__
    void hessian(const double* x, double* H) const {
        static_cast<const Derived*>(this)->hessian_impl(x, H);
    }
};

} // namespace numa::cuda
