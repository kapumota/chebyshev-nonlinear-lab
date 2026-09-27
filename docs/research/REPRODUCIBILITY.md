### Reproducibilidad

#### Dependencias de referencia

Entorno Linux recomendado:

```bash
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    ninja-build \
    git \
    zip \
    libomp-dev \
    openmpi-bin \
    libopenmpi-dev
```

#### CPU y OpenMP

```bash
cmake \
    -S . \
    -B build \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DNUMA_BUILD_BENCH=OFF

cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Resultado de referencia de R0:

```text
21/21 tests correctos
```

#### MPI

```bash
cmake \
    -S . \
    -B build-mpi \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DNUMA_BUILD_MPI=ON \
    -DNUMA_BUILD_BENCH=OFF

cmake --build build-mpi -j"$(nproc)"
ctest --test-dir build-mpi -L mpi --output-on-failure
```

CTest registra configuraciones con 1, 2, 3 y 8 procesos.

#### Sanitizers

```bash
./scripts/run_sanitizers.sh
```

#### CUDA

La ausencia de GPU no debe transformarse en un `PASS` ni en un `FAIL` ficticio.

Registrar:

```text
NOT_EXECUTED_NO_GPU
```

cuando corresponda.

#### Captura del entorno

Antes de una campaña cientifica registrar:

```bash
uname -a
lscpu
free -h
g++ --version
cmake --version
mpicxx --version
mpirun --version
```

Para GPU, cuando exista:

```bash
nvidia-smi
nvcc --version
```

#### Reproduccion limpia de R0

El freeze `baseline-v0.1.0` requiere una prueba desde otro directorio:

```bash
rm -rf /tmp/chebyshev-r0-reproduction

git clone \
    /ruta/al/repositorio/chebyshev-nonlinear-lab \
    /tmp/chebyshev-r0-reproduction

cd /tmp/chebyshev-r0-reproduction
git switch r0/reconciliacion-baseline
```

Despues se repiten los procedimientos CPU/OpenMP y MPI anteriores.

La reproduccion debe partir de un arbol limpio y no depender de `build/`, caches ni archivos no versionados del workspace original.
