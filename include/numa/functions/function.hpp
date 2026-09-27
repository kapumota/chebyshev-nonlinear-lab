#pragma once
#include <array>
#include "../core/vector.hpp"
#include "../core/matrix.hpp"

namespace numa {

template <std::size_t Dim>
struct IFunction {
    virtual ~IFunction() = default;
    virtual Vec<Dim> eval(const Vec<Dim>& x) const = 0;
    virtual Matrix<Dim, Dim> jacobian(const Vec<Dim>& x) const = 0;
    virtual Matrix<Dim, Dim * Dim> hessian(const Vec<Dim>& x) const = 0;
    virtual double lipschitz_k() const = 0;

    // Acceso por fila, para distribucion MPI (Proyecto 3, ver
    // include/numa/mpi/chebyshev_mpi.hpp).
    //
    // Implementacion por defecto: llama a eval()/jacobian()/hessian()
    // completos y extrae la fila i. No ahorra trabajo, pero mantiene
    // compatibilidad con cualquier IFunction existente sin tocarlo.
    //
    // Los problemas donde el costo esta en la evaluacion (Chandrasekhar,
    // o el reactor generalizado del Proyecto 7) deben sobreescribir
    // estos tres metodos con una version que calcule solo la fila
    // pedida, en O(N) en vez de O(N^2)/O(N^3). Ver el override real en
    // include/numa/problems/chandrasekhar.hpp de este mismo paquete.
    virtual double row(std::size_t i, const Vec<Dim>& x) const {
        return eval(x)[i];
    }
    virtual std::array<double, Dim> jacobian_row(std::size_t i, const Vec<Dim>& x) const {
        return jacobian(x)[i];
    }
    virtual std::array<double, Dim * Dim> hessian_row(std::size_t i, const Vec<Dim>& x) const {
        return hessian(x)[i];
    }
};

} // namespace numa
