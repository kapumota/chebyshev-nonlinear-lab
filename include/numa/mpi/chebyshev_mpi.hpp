#pragma once
#ifdef NUMA_HAS_MPI
#include <mpi.h>
#endif
#include <algorithm>
#include <vector>
#include "../core/vector.hpp"
#include "../core/matrix.hpp"
#include "../functions/function.hpp"
#include "../linalg/gauss_seidel.hpp"
#include "../linalg/tri_sup.hpp"

// Chebyshev distribuido con MPI.
//
// CORRECCION respecto a la version original (stub del Proyecto 15):
// el stub calculaba F, J y H completos, identicos, en cada rank
// (llamaba a f.eval(x) / f.jacobian(x) sobre el vector x completo,
// no sobre un bloque), y luego ensamblaba con
// MPI_Allreduce(..., MPI_SUM, ...). Como todos los ranks aportaban
// el mismo valor, el Allreduce multiplicaba el resultado por `size`
// en vez de devolver el valor real. El metodo no distribuia trabajo
// y ademas era numericamente incorrecto para size > 1.
//
// Esta version reparte las FILAS de F, J y H entre ranks. Cada rank
// calcula solo su bloque via IFunction::row/jacobian_row/hessian_row
// (ver include/numa/functions/function.hpp), y el resultado se
// ensambla con MPI_Allgatherv, que concatena bloques distintos en
// vez de sumarlos. El paso lineal (gauss_seidel + tri_sup) se
// resuelve de forma redundante en cada rank sobre la matriz ya
// ensamblada: la ganancia no esta en el algebra lineal (barata para
// N moderado), sino en evitar que cada rank recalcule F/J/H
// completos cuando esas evaluaciones son el costo dominante.
//
// Requiere que IFunction<N> tenga row()/jacobian_row()/hessian_row()
// sobreescritos con costo real por fila (ver el override en
// include/numa/problems/chandrasekhar.hpp de este mismo paquete).
// Con la implementacion por defecto de la interfaz (que llama a
// eval()/jacobian()/hessian() completos y solo extrae la fila) este
// archivo sigue siendo correcto, pero no gana nada en rendimiento
// frente a la version secuencial: el ahorro depende enteramente de
// que el problema tenga overrides por fila reales.

namespace numa::mpi {

template <std::size_t N>
struct MpiResult {
    Vec<N> x;
    std::size_t iterations;
    double er_abs;
    int rank, size;
};

namespace detail {

// Reparte N filas entre `size` ranks lo mas parejo posible. Los
// primeros (N % size) ranks reciben una fila extra.
inline void row_range(int N_total, int rank, int size,
                       int& row_start, int& local_count)
{
    const int base = N_total / size;
    const int rem  = N_total % size;
    row_start   = rank * base + std::min(rank, rem);
    local_count = base + (rank < rem ? 1 : 0);
}

} // namespace detail

template <std::size_t N>
MpiResult<N> chebyshev_mpi(const IFunction<N>& f,
                            Vec<N> x,
                            int rank, int size,
                            double tol = 1e-12,
                            std::size_t max_iter = 20)
{
    int row_start = 0, local_count = static_cast<int>(N);
    detail::row_range(static_cast<int>(N), rank, size, row_start, local_count);
    const int row_end = row_start + local_count;

#ifdef NUMA_HAS_MPI
    // Tablas de counts/displs para MPI_Allgatherv, precalculadas una
    // vez fuera del bucle de iteracion. counts_h/displs_h se dejan
    // listos para el caso en que se distribuya tambien el ensamblado
    // del hessiano completo; en esta version el hessiano se reduce a
    // un vector de tamano N (ver mas abajo), asi que solo hacen
    // falta counts/displs (para F, w) y counts_j/displs_j (para J,
    // que tiene N columnas por fila).
    std::vector<int> counts(size), displs(size), counts_j(size), displs_j(size);
    for (int r = 0; r < size; ++r) {
        int rs = 0, lc = static_cast<int>(N);
        detail::row_range(static_cast<int>(N), r, size, rs, lc);
        counts[r]   = lc;
        displs[r]   = rs;
        counts_j[r] = lc * static_cast<int>(N);
        displs_j[r] = rs * static_cast<int>(N);
    }
#endif

    double er = 0.0;
    std::size_t k = 0;
    Vec<N> x_ref = x;

    while (k == 0 || (er >= tol && k <= max_iter)) {
        // --- Ensamblado distribuido de F y J ---
        Vec<N> F{};
        Matrix<N, N> J{};
        for (int i = row_start; i < row_end; ++i) {
            F[i] = f.row(static_cast<std::size_t>(i), x);
            J[i] = f.jacobian_row(static_cast<std::size_t>(i), x);
        }

#ifdef NUMA_HAS_MPI
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, F.data(),
                        counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, J.data()->data(),
                        counts_j.data(), displs_j.data(), MPI_DOUBLE, MPI_COMM_WORLD);
#endif

        Matrix<N, N> J1 = J; // copia: gauss_seidel destruye la matriz de entrada
        Vec<N> F1 = F;
        gauss_seidel(J1, F1);
        Vec<N> z = tri_sup(J1, F1);

        // --- Termino cuadratico de Chebyshev, tambien distribuido ---
        Vec<N> w{};
        for (int i = row_start; i < row_end; ++i) {
            auto Hi = f.hessian_row(static_cast<std::size_t>(i), x);
            double s = 0.0;
            for (std::size_t a = 0; a < N; ++a)
                for (std::size_t b = 0; b < N; ++b)
                    s += Hi[a * N + b] * z[a] * z[b];
            w[i] = s;
        }
#ifdef NUMA_HAS_MPI
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, w.data(),
                        counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
#endif

        // Segunda resolucion lineal, con una copia fresca de J (la
        // primera copia J1 ya quedo triangularizada por gauss_seidel
        // y no sirve para un segundo sistema). Este es el paso que
        // faltaba en el borrador anterior: gauss_seidel(J2, w) tiene
        // que ejecutarse ANTES de tri_sup(J2, w), porque tri_sup
        // asume una matriz ya triangular superior.
        Matrix<N, N> J2 = J;
        Vec<N> w2 = w;
        gauss_seidel(J2, w2);
        Vec<N> u = tri_sup(J2, w2);

        Vec<N> dx{};
        for (std::size_t i = 0; i < N; ++i) dx[i] = z[i] + 0.5 * u[i];

        x = x - dx;
        er = norm_inf(x - x_ref);
        x_ref = x;
        ++k;
    }

    return { x, k, er, rank, size };
}

} // namespace numa::mpi
