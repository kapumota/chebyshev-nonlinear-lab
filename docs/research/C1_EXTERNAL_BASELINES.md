### Baselines externos obligatorios de C1

#### Estado

```text
REQUIRED_NOT_YET_FROZEN
```

El Paper 1 no se cierra únicamente con implementaciones internas. Antes del freeze confirmatorio debe existir al menos un baseline externo mantenido por terceros, ejecutado con versión y opciones exactas registradas.

#### Baseline requerido

El baseline externo principal es PETSc `SNESNEWTONLS`, Newton con line search. PETSc usa `KSP` para las resoluciones lineales internas. La versión, el tipo de `KSP`, el precondicionador, tolerancias, line search y todas las opciones efectivas deben quedar congeladas antes de C1.

No se declara que una configuración PETSc sea equivalente a nuestro Newton-Krylov. Su función es impedir que el resultado dependa únicamente de una implementación interna potencialmente favorable o desfavorable.

#### Requisitos mínimos del manifest

Cada ejecución externa debe quedar asociada a un manifest JSON con:

```text
schema_version
solver
solver_version
nonlinear_method
linear_method
preconditioner
matrix_mode
options
repository_commit
executable_sha256
build
hardware
command
problems
created_utc
```

El manifest no contiene resultados científicos resumidos. Los resultados se guardan aparte y se relacionan mediante hash y configuración.

#### Reglas de equidad

1. El baseline externo recibe las mismas instancias científicas y condiciones iniciales que el método interno cuando la representación lo permita.
2. No se eliminan instancias por fallo de PETSc ni por fallo de un método interno.
3. Las tolerancias no lineales deben corresponder al mismo criterio de éxito preregistrado.
4. El tuning del baseline externo se realiza sólo en calibración.
5. La configuración seleccionada se congela antes de confirmación.
6. El costo principal sigue siendo `W_eq`, tasa de éxito y `t_total`; las métricas nativas de PETSc se preservan como diagnóstico.
7. La evidencia confirmatoria no puede usarse para cambiar el `KSP`, `PC`, line search o tolerancias.

#### Gate

`READY_FOR_C1` requiere un manifest válido y versionado. La sola presencia de PETSc en el sistema no cierra el gate. También deben existir scripts o comandos reproducibles que generen los resultados externos desde las instancias congeladas.
