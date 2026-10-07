### Chebyshev Nonlinear Lab

Framework de investigación reproducible en C++20 para métodos iterativos no lineales, con énfasis en métodos de orden superior, computación dispersa, ejecución paralela y evaluación a gran escala.

El repositorio contiene un baseline funcional con Newton, Newton por diferencias finitas, Broyden, Newton-Kantorovich y Chebyshev, junto con rutas OpenMP, MPI y CUDA. Los problemas pequeños actualmente incluidos se conservan como pruebas numéricas, golden tests y regresión. No constituyen la escala experimental objetivo de los futuros trabajos científicos.

#### Estado del baseline

El baseline reconciliado dispone de:

- C++20 y CMake
- tests con Catch2
- OpenMP
- verificación MPI automatizada con CMake y CTest para 1, 2, 3 y 8 procesos
- Google Benchmark
- ruta CUDA preservada para validación posterior con hardware adecuado
- scripts de compilación, tests, sanitizers y empaquetado

La línea de investigación posterior al baseline se orienta a sistemas no lineales dispersos de gran dimensión, álgebra lineal iterativa, operaciones matrix-free, derivadas direccionales de segundo orden, estrategias adaptativas y evaluación reproducible sobre benchmarks escalables y aplicaciones públicas.

#### Alcance de los casos pequeños

Los casos 2x2 y 8x8 cumplen exclusivamente funciones de validación:

- comprobar fórmulas y criterios de convergencia
- detectar regresiones
- contrastar implementaciones secuenciales y paralelas
- verificar reproducibilidad entre configuraciones

La evaluación científica futura se regirá por `docs/research/LARGE_SCALE_PROTOCOL.md`.

#### Compilación CPU y OpenMP

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

#### Compilación y tests MPI

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

También puede utilizarse:

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

Los benchmarks actuales pertenecen al baseline. Los protocolos de gran escala se incorporarán en fases posteriores.

#### CUDA

La ruta CUDA se conserva en `include/numa/cuda/` y `apps/demo_chebyshev_gpu.cu`. Su validación experimental no forma parte del freeze inicial cuando no existe una GPU NVIDIA disponible.

Un estado no ejecutado por ausencia de hardware debe registrarse como:

```text
NOT_EXECUTED_NO_GPU
```

#### Documentación de investigación

La documentación principal se encuentra en `docs/research/`:

- `BASELINE.md`, estado congelado del software
- `PROVENANCE.md`, procedencia y trazabilidad del baseline
- `ROADMAP.md`, transición hacia investigación de gran escala
- `REPRODUCIBILITY.md`, procedimiento de reproducción
- `LARGE_SCALE_PROTOCOL.md`, escalas y métricas experimentales
- `STATE_OF_THE_ART.md`, gate de auditoría bibliográfica y novedad
- `STYLE_GUIDE.md`, convenciones del repositorio
- `PUBLIC_RELEASE_CHECKLIST.md`, comprobaciones previas a hacer público el repositorio

#### Líneas de investigación

Las carpetas `experiments/`, `manuscript/` y `artifacts/` reservan espacio para varias líneas de investigación. Los títulos y el número final de papers no se consideran congelados. Antes de desarrollar R1 se realizará una auditoría adversarial del estado del arte que podrá producir decisiones `GO`, `REFRAME`, `MERGE` o `KILL` para cada línea.

#### Reproducibilidad

El baseline original permanece congelado mediante el tag `source-v0.0.0`. El siguiente freeze estable será `baseline-v0.1.0` una vez completada una reproducción desde clon limpio.

Consulte `docs/research/REPRODUCIBILITY.md` para el protocolo completo.

#### Licencia

El código propio del repositorio se distribuye bajo la licencia MIT, salvo componentes externos que mantienen sus respectivas licencias.
