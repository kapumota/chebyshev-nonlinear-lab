// Verifica automáticamente que Chebyshev MPI conserva la solución de
// referencia para Chandrasekhar<8> con cualquier número de procesos usado
// por CTest.
#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdio>

#include <numa/mpi/chebyshev_mpi.hpp>
#include <numa/problems/chandrasekhar.hpp>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    using namespace numa;

    numa::problems::Chandrasekhar<8> problem;
    Vec<8> x0{};
    for (auto& value : x0) {
        value = 1.0;
    }

    const auto result = numa::mpi::chebyshev_mpi<8>(
        problem, x0, rank, size, 1e-12, 20);

    constexpr std::array<double, 8> reference{
        1.021719731462,
        1.073186381734,
        1.125724893657,
        1.169753312169,
        1.203071751305,
        1.226490874633,
        1.241524600594,
        1.249448516693
    };

    constexpr double solution_tolerance = 1e-10;
    constexpr double residual_tolerance = 1e-10;

    bool local_ok = true;
    for (std::size_t i = 0; i < reference.size(); ++i) {
        if (std::abs(result.x[i] - reference[i]) > solution_tolerance) {
            local_ok = false;
        }
    }

    const auto residual = problem.eval(result.x);
    if (norm_inf(residual) > residual_tolerance) {
        local_ok = false;
    }

    int local_status = local_ok ? 1 : 0;
    int global_status = 0;
    MPI_Allreduce(
        &local_status,
        &global_status,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD);

    if (rank == 0) {
        std::printf(
            "MPI procesos=%d, iteraciones=%zu, error_paso=%.6e, residual=%.6e\n",
            size,
            result.iterations,
            result.er_abs,
            norm_inf(residual));

        for (std::size_t i = 0; i < reference.size(); ++i) {
            std::printf("x[%zu]=%.12f\n", i, result.x[i]);
        }

        std::printf(
            "%s\n",
            global_status == 1
                ? "[ok] La solución MPI coincide con la referencia"
                : "[FALLA] La solución MPI no coincide con la referencia");
    }

    MPI_Finalize();
    return global_status == 1 ? 0 : 1;
}
