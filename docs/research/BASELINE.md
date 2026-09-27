### Baseline reproducible

#### Identidad

Repositorio: `chebyshev-nonlinear-lab`

Estado de trabajo: R0, reconciliacion y freeze del baseline.

Tag de origen inmutable:

```text
source-v0.0.0
```

#### Estado funcional

El baseline reconciliado incluye:

- Newton
- Newton con diferencias finitas
- Broyden
- Newton-Kantorovich
- Chebyshev
- variantes OpenMP
- verificacion MPI integrada en CMake y CTest
- ruta CUDA preservada
- tests unitarios, golden y de regresion

Los sistemas 2x2 y 8x8 se clasifican como casos de validacion. No constituyen evidencia experimental de gran escala.

#### Gates confirmados

```text
Build C++20             PASS
OpenMP                   PASS
CTest CPU/OpenMP         21/21 PASS
Verificacion secuencial  PASS
MPI manual 1,2,3,8       PASS
MPI integrado en CMake   IMPLEMENTED
CUDA                     NOT_EXECUTED_NO_GPU
```

El gate MPI integrado se considera cerrado cuando CTest reporte 4/4 tests MPI correctos desde `build-mpi`.

#### Correcciones reconciliadas

El baseline importado presentaba defectos detectados por los tests. La rama R0 los corrige sin alterar el tag `source-v0.0.0`.

Commits de reconciliacion actualmente registrados:

```text
7b3f753  Corregir recurrencias y criterios de convergencia del baseline
bc5fa60  Integrar MPI en CMake y automatizar su verificacion
```

#### Freeze siguiente

El tag:

```text
baseline-v0.1.0
```

solo se crea despues de una reproduccion completa desde un clon limpio y una auditoria `PUBLIC-READY`.
