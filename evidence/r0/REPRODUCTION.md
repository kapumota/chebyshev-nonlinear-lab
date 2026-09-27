### Reproducción limpia del baseline

#### Estado

R0-G10: PASS.

La reproducción se realizó desde un clon limpio de la rama:

`r0/reconciliacion-baseline`

Candidato reproducido:

`cef2e57`

#### CPU y OpenMP

La configuración y compilación se realizaron desde cero con CMake y Ninja.

Resultado:

`21/21 PASS`

#### MPI

Se realizó una configuración independiente con:

`NUMA_BUILD_MPI=ON`

Se verificaron:

- `mpi_np1`
- `mpi_np2`
- `mpi_np3`
- `mpi_np8`

Resultado:

`4/4 PASS`

#### CUDA

Estado:

`NOT_EXECUTED_NO_GPU`

La ausencia de ejecución CUDA no se interpreta como PASS ni como FAIL.

#### Packaging

`scripts/package.sh` se ejecutó correctamente desde el clon limpio.

SHA-256 del paquete generado durante esta reproducción:

`4b4e53620ad2bb49be75eb701c453bff45270c8ae76ed29bacefc81cf8281d50`

Este checksum corresponde al candidato reproducido antes del merge final y antes del tag `baseline-v0.1.0`.

El paquete definitivo deberá regenerarse después del merge y de la creación del tag.

#### Conclusión

El baseline puede reconstruirse desde un clon limpio y reproduce correctamente las rutas CPU, OpenMP y MPI verificadas durante R0.
