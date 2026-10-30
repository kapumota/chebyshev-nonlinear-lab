### Runner PETSc real para C1

#### Alcance de 0020b-v1

Esta entrega integra un ejecutable PETSc real para B1 Bratu. No amplía todavía la campaña a B3a, B5 o Bext y no constituye el freeze del baseline externo.

El sistema no lineal es el mismo objeto `Bratu3D` usado por `run_campaign`. PETSc recibe:

```text
F(x)        -> Bratu3D::eval
J(x) v      -> Bratu3D::jvp mediante MatShell
Jacobi      -> Bratu3D::jacobi_diag mediante PCSHELL
```

No se materializa un Jacobiano y no se usa `-snes_mf` como camino principal. Por tanto, la acción del Jacobiano corresponde al Jv analítico matrix-free del repositorio y no a diferencias finitas internas de PETSc.

#### Configuración efectiva

El ejecutable fija valores iniciales reproducibles para:

```text
SNESNEWTONLS
SNESLINESEARCHBT
KSPGMRES
MatShell
PCSHELL o PCNONE según --precond
```

Luego aplica `SNESSetFromOptions`. Las opciones PETSc pueden modificar SNES, line search, KSP y tolerancias. El operador `MatShell` y el precondicionador solicitado por `--precond` se reinstalan después para impedir que una opción externa sustituya silenciosamente F, Jv o la política de precondicionamiento de la instancia.

El ejecutable consulta después los tipos y tolerancias efectivos y los emite entre `PETSC_BASELINE_BEGIN` y `PETSC_BASELINE_END`. El manifest se construye a partir de esos valores observados, no de una plantilla rellenada manualmente.

#### Contabilidad

El CSV utiliza el mismo esquema crudo de C1. Para PETSc se registran:

```text
n_F          llamadas reales al callback F
n_JVP        MatMult reales del MatShell
n_prec       aplicaciones reales del PCSHELL
n_D2         0
n_lin        SNESGetLinearSolveIterations
t_setup      actualización de la diagonal Jacobi
t_F
t_JVP
t_prec
t_total      tiempo de SNESSolve
success      SNESConvergedReason > 0
```

`n_JVP` sigue siendo una métrica mecanística. No sustituye `W_eq`, tasa de éxito ni `t_total`.

#### Captura de procedencia

`scripts/capture_petsc_baseline.py` exige árboles Git limpios para el proyecto y PETSc, ejecuta el baseline con un hilo, calcula el SHA-256 del ejecutable y obtiene la configuración efectiva emitida por PETSc. También registra commit y tag de PETSc, comando, compilación, CPU y hashes de la salida cruda.

El manifest producido por 0020b-v1 sólo declara `B1`. No debe promoverse a manifest final de C1 hasta que el mismo runner cubra B3a, B5 y Bext, se ejecute la calibración externa y se congele la configuración seleccionada.
