### Chebyshev Nonlinear Lab

Framework de investigacion reproducible en C++20 para metodos iterativos no lineales, con enfasis en metodos de orden superior, computacion dispersa, ejecucion paralela y evaluacion a gran escala.

El repositorio contiene un baseline funcional con Newton, Newton por diferencias finitas, Broyden, Newton-Kantorovich y Chebyshev, junto con rutas OpenMP, MPI y CUDA. Los problemas pequenos actualmente incluidos se conservan como pruebas numericas, golden tests y regresion. No constituyen la escala experimental objetivo de los futuros trabajos cientificos.

#### Estado del baseline

El baseline reconciliado dispone de:

- C++20 y CMake
- tests con Catch2
- OpenMP
- verificacion MPI automatizada con CMake y CTest para 1, 2, 3 y 8 procesos
- Google Benchmark
- ruta CUDA preservada para validacion posterior con hardware adecuado
- scripts de compilacion, tests, sanitizers y empaquetado

La linea de investigacion posterior al baseline se orienta a sistemas no lineales dispersos de gran dimension, algebra lineal iterativa, operaciones matrix-free, derivadas direccionales de segundo orden, estrategias adaptativas y evaluacion reproducible sobre benchmarks escalables y aplicaciones publicas.

#### Alcance de los casos pequenos

Los casos 2x2 y 8x8 cumplen exclusivamente funciones de validacion:

- comprobar formulas y criterios de convergencia
- detectar regresiones
- contrastar implementaciones secuenciales y paralelas
- verificar reproducibilidad entre configuraciones

La evaluacion cientifica futura se regira por `docs/research/LARGE_SCALE_PROTOCOL.md`.

#### Compilacion CPU y OpenMP

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

#### Compilacion y tests MPI

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

Tambien puede utilizarse:

```bash
./scripts/build_mpi.sh
```

#### Sanitizers

```bash
./scripts/run_sanitizers.sh
```

#### Benchmarks

```bash
./scripts/run_bench.sh
```

Los benchmarks actuales pertenecen al baseline. Los protocolos de gran escala se incorporaran en fases posteriores.

#### CUDA

La ruta CUDA se conserva en `include/numa/cuda/` y `apps/demo_chebyshev_gpu.cu`. Su validacion experimental no forma parte del freeze inicial cuando no existe una GPU NVIDIA disponible.

Un estado no ejecutado por ausencia de hardware debe registrarse como:

```text
NOT_EXECUTED_NO_GPU
```

#### Documentacion de investigacion

La documentacion principal se encuentra en `docs/research/`:

- `BASELINE.md`, estado congelado del software
- `PROVENANCE.md`, procedencia y trazabilidad del baseline
- `ROADMAP.md`, transicion hacia investigacion de gran escala
- `REPRODUCIBILITY.md`, procedimiento de reproduccion
- `LARGE_SCALE_PROTOCOL.md`, escalas y metricas experimentales
- `STATE_OF_THE_ART.md`, gate de auditoria bibliografica y novedad
- `STYLE_GUIDE.md`, convenciones del repositorio
- `PUBLIC_RELEASE_CHECKLIST.md`, comprobaciones previas a hacer publico el repositorio

#### Lineas de investigacion

Las carpetas `experiments/`, `manuscript/` y `artifacts/` reservan espacio para varias lineas de investigacion. Los titulos y el numero final de papers no se consideran congelados. Antes de desarrollar R1 se realizara una auditoria adversarial del estado del arte que podra producir decisiones `GO`, `REFRAME`, `MERGE` o `KILL` para cada linea.

#### Reproducibilidad

El baseline original permanece congelado mediante el tag `source-v0.0.0`. El siguiente freeze estable sera `baseline-v0.1.0` una vez completada una reproduccion desde clon limpio.

Consulte `docs/research/REPRODUCIBILITY.md` para el protocolo completo.

#### Licencia

El codigo propio del repositorio se distribuye bajo la licencia MIT, salvo componentes externos que mantienen sus respectivas licencias.
