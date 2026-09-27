### chebyshev-solver

Implementación en **C++20** de los 5 algoritmos del Anexo I de la tesis
*"Solución de un sistema no lineal algebraico por optimización numérica"*
(Leopoldo Paredes Soria, UNI, 2021).

Incluye versiones **secuenciales**, **paralelas (OpenMP)** y **GPU (CUDA)**,
más un stub **MPI** para sistemas de mayor dimensión.

#### Algoritmos

| # | Algoritmo | Archivo | Caso |
|---|---|---|---|
| 1 | Newton clásico | `newton.hpp` | 2×2 |
| 2 | Newton con diferencias finitas | `newton_df.hpp` | 2×2 |
| 3 | Broyden | `broyden.hpp` | 2×2 |
| 4 | Newton-Kantorovich | `kantorovich.hpp` | 8×8 |
| 5 | Chebyshev iterativo | `chebyshev.hpp` | 8×8 |

**Paralelos:** `newton_omp.hpp`, `chebyshev_omp.hpp` (OpenMP)
**GPU:** `cuda/chebyshev_kernel.cuh` (CUDA)
**Distribuido:** `mpi/chebyshev_mpi.hpp` (MPI, stub)

#### Requisitos

- CMake ≥ 3.20
- Compilador C++20 (GCC ≥ 11, Clang ≥ 14)
- OpenMP (opcional, recomendado)
- CUDA Toolkit ≥ 11 (opcional)
- MPI (opcional)

#### Build

```bash
# Build secuencial + OpenMP
./scripts/build.sh

# Build con CUDA
./scripts/build_cuda.sh

# Build con MPI
./scripts/build_mpi.sh

# Build con sanitizers (ASan + UBSan)
./scripts/run_sanitizers.sh
# ============================================================
# .gitignore
# ============================================================
cat > .gitignore << 'NUMA_EOF'
build/
build-*/
*.o
*.a
*.so
*.zip
*.data
.cache/
compile_commands.json
.cmake/
