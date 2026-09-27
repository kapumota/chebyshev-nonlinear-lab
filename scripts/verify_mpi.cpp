// Verifica que chebyshev_mpi da el mismo resultado sin importar cuantos
// ranks se usen (prueba de correccion del bug de Allreduce/tri_sup).
//
// Compilar: mpicxx -std=c++20 -Iinclude -DNUMA_HAS_MPI -O2 scripts/verify_mpi.cpp -o /tmp/verify_mpi
// Correr:   mpirun --oversubscribe -np 3 /tmp/verify_mpi
#include <mpi.h>
#include <cstdio>
#include <numa/problems/chandrasekhar.hpp>
#include <numa/mpi/chebyshev_mpi.hpp>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    using namespace numa;
    numa::problems::Chandrasekhar<8> ch;
    Vec<8> x0{}; for (auto& v : x0) v = 1.0;

    auto r = numa::mpi::chebyshev_mpi<8>(ch, x0, rank, size, 1e-12, 20);

    if (rank == 0) {
        std::printf("size=%d  iter=%zu  er=%e\n", size, r.iterations, r.er_abs);
        for (int i = 0; i < 8; ++i) std::printf("x[%d]=%.12f\n", i, r.x[i]);
        std::printf("\nCorre esto mismo con -np 1, -np 2, -np 3, -np 8 y compara:\n");
        std::printf("todas las filas de x deben coincidir en al menos 10 decimales.\n");
    }
    MPI_Finalize();
    return 0;
}
