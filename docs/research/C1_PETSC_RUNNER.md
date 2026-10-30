### Runner PETSc externo para C1

`petsc_baseline` ejecuta PETSc sobre los mismos operadores del repositorio y
no sobre discretizaciones reimplementadas. Las cuatro familias admitidas son:

- B1: Bratu 3D, `MatShell` sobre `Bratu3D::jvp` y Jacobi mediante `PCSHELL`.
- B5: Chandrasekhar dinámico, `MatShell` sobre `ChandrasekharOp::jvp`, sin PC.
- B3a: Burgers 2D manufacturado de control, `MatShell` y MG mediante `PCSHELL`.
- Bext: Brusselator 1D no manufacturado, `MatShell` y Jacobi por bloques mediante `PCSHELL`.

El solver externo se fija a `SNESNEWTONLS`, line search `bt` y `KSPGMRES`.
Después de `SNESSetFromOptions`, el ejecutable consulta la configuración
efectiva y aborta si esas tres decisiones fueron sustituidas. El tipo de PC
queda fijado por la familia y por `--precond`.

### Starts materializados

`natural` usa exactamente el estado natural de cada problema. Para
`perturbed`, `continuation` y `coarse`, PETSc no genera un vector propio: exige
`--start-file` y consume los mismos valores materializados que la campaña
interna. De este modo la comparación no depende de que dos generadores
pseudoaleatorios produzcan accidentalmente el mismo vector.

Los campos `start_id`, `distance`, `rho_x` y `seed` se copian al CSV y forman
parte de la procedencia. La construcción de esos starts pertenece al protocolo
C1, no al baseline PETSc.

### Instrumentación

Cada fila conserva `n_F`, `n_JVP`, `n_prec`, iteraciones no lineales y
lineales, fallos lineales, tiempos de primitivas, `t_total`, residual relativo,
RSS, commit y configuración efectiva. PETSc tiene `n_D2=0`.

`n_JVP` es una variable mecanicista y no un objetivo de selección. La
calibración posterior mantiene la regla global `success -> W_eq -> t_total`.

### Alcance

B3a sigue siendo un benchmark controlado manufacturado. Bext proporciona el
caso público no manufacturado. Ningún resultado de este runner autoriza a
retirar regiones donde Newton-Krylov gane ni a modificar la malla de
confirmación.
