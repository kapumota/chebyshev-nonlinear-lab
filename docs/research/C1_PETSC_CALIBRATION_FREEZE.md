### Calibración y freeze de PETSc para C1

#### Estado

Este protocolo cierra el paso entre el runner PETSc de cuatro familias y el gate `READY_FOR_C1`. No contiene resultados confirmatorios y no autoriza ejecutar C1 antes de cerrar `0021-v4`, `0022-v4` y `0023-v4`.

#### Invariantes del mecanismo

Las cuatro familias usan el mismo mecanismo externo:

```text
SNESNEWTONLS
SNESLINESEARCHBT
KSPGMRES
MatShell analytic Jv
un proceso
un hilo
```

El operador no se vuelve a implementar dentro de PETSc. El callback `F`, el `Jv` analítico y el precondicionador de aplicación proceden de los mismos operadores del repositorio ya validados por `0020c-v1`.

El precondicionador puede formar parte del espacio de calibración únicamente entre variantes ya implementadas por la familia. No se admite introducir un nuevo precondicionador después de observar C1.

#### Separación estricta de fases

Los smoke tests previos permanecen `exploratory`. No entran en el ranking de PETSc.

La calibración PETSc:

```text
usa exclusivamente phase=calibration
usa parámetros y starts reservados para calibración
usa un commit de código sin cambios rastreados
usa un plan JSON fijado antes de ejecutar
conserva todos los candidatos, incluidos los perdedores
no usa ninguna fila confirmatoria
```

`run_petsc_calibration.py` rechaza opciones reservadas dentro de las instancias y fuerza la fase, los límites y el núcleo PETSc desde el plan.

#### Presupuesto y espacio de candidatos

El plan declara `expected_candidate_count`. Ese valor debe ser el mismo presupuesto máximo de configuraciones fijado para los métodos comparados. Cada familia PETSc recibe exactamente ese número de candidatos.

Los límites de exactitud y trabajo son presupuesto, no hiperparámetros competitivos:

```text
rtol
max_nl
maxlin
warmup
reps
threads
```

El núcleo también permanece fijo:

```text
snes_type          newtonls
linesearch_type    bt
ksp_type           gmres
matrix_mode        MatShell analytic Jv
```

La calibración puede decidir, según el plan previo:

```text
precond            sólo entre variantes ya implementadas para la familia
ksp_rtol           tolerancia lineal PETSc candidata
gmres_restart      restart candidato
```

No se añaden dimensiones de tuning después de observar qué configuración gana.

#### Esquema del plan

Las rutas se interpretan desde la raíz del repositorio. El siguiente fragmento sólo documenta el esquema y es deliberadamente no ejecutable. No fija el presupuesto real de candidatos ni sustituye las mallas de calibración ya acordadas:

```json
{
  "schema_version": 1,
  "campaign": "C1",
  "phase": "calibration",
  "source_commit": "SHA_COMPLETO_DEL_COMMIT_DE_CODIGO",
  "executable": "build-petsc/petsc_baseline",
  "expected_candidate_count": 4,
  "fixed": {
    "snes_type": "newtonls",
    "linesearch_type": "bt",
    "ksp_type": "gmres",
    "matrix_mode": "MatShell analytic Jv",
    "rtol": 1e-8,
    "max_nl": 50,
    "maxlin": 200,
    "warmup": 1,
    "reps": 7,
    "threads": 1
  },
  "families": {
    "B1": {
      "instances": [
        {
          "id": "identificador_calibracion_B1",
          "args": [
            "--problem", "bratu", "--family", "B1",
            "--n", "31", "--param", "VALOR_DE_CALIBRACION",
            "--start", "natural"
          ]
        }
      ],
      "candidates": [
        {
          "id": "p01",
          "precond": "jacobi",
          "ksp_rtol": 1e-5,
          "restart": 30
        }
      ]
    },
    "B5": {"instances": ["..."], "candidates": ["..."]},
    "B3a": {"instances": ["..."], "candidates": ["..."]},
    "Bext": {"instances": ["..."], "candidates": ["..."]}
  }
}
```

El fragmento no es JSON ejecutable. El archivo real debe materializar el presupuesto y las mallas ya acordadas antes de ejecutarse. El `source_commit` identifica el código y puede coexistir con un plan todavía no rastreado, porque el runner rechaza cambios de código rastreados pero permite artefactos nuevos cuyo SHA-256 queda capturado antes de la primera corrida.

#### Flujo obligatorio

```text
plan PETSc previo
        |
        v
run_petsc_calibration.py
        |
        v
PETSc raw calibration CSV + run manifest
        |
        +------------------------------+
        |                              |
        v                              v
freeze_primitive_costs.py         auditoría de costos
        |                              |
        v                              |
primitive_costs frozen <---------------+
        |
        v
compute_weq.py sobre PETSc raw
        |
        v
PETSc calibration CSV + W_eq
        |
        v
freeze_petsc_calibration.py
        |
        +--> ranking completo
        +--> petsc_config frozen
        +--> freeze manifest con SHA-256
        |
        v
verify_petsc_freeze.py
```

La tabla de costos conserva la definición ya congelada:

```text
W_eq =
    (n_F*c_F
     + n_JVP*c_JVP
     + n_D2*c_D2
     + n_prec*c_prec
     + t_setup)
    / c_JVP
```

PETSc tiene `n_D2=0`, pero no recibe una definición especial de `W_eq`.

#### Regla de selección

Se aplica exactamente la misma semántica de tuning ya fijada para C1:

```text
1  maximizar tasa de éxito
2  minimizar mediana de W_eq entre instancias exitosas
3  minimizar mediana de t_total entre instancias exitosas
```

Las repeticiones temporales se condensan por mediana antes de comparar configuraciones. Una configuración fallida en una instancia no recibe un costo finito inventado para esa instancia. La tasa de éxito se evalúa antes que el costo.

Si dos candidatos empatan exactamente en los tres criterios, el identificador del candidato sólo resuelve el orden de manera determinista. No constituye un cuarto objetivo científico.

#### Freeze

El freeze produce:

```text
ranking PETSc completo
configuración seleccionada por B1, B5, B3a y Bext
hash del plan previo
hash del CSV crudo
hash del CSV con W_eq
hash de la tabla de costos
hash de archivos de procedencia adicionales
hash del ranking
hash de la configuración congelada
hash del selector
source_commit de calibración
```

El estado válido es:

```text
FROZEN_CALIBRATION_ONLY
```

No significa `READY_FOR_C1`. Ese estado sólo puede ser concedido por `0023-v4` después de cerrar el análisis estadístico, el preregistro final, los starts, TF2018, hardware, costos y manifests.

#### Prohibiciones después del freeze

Después del freeze PETSc no se permite:

```text
cambiar el candidato seleccionado porque otro rinde mejor en C1
cambiar precondicionador, ksp_rtol o restart tras mirar confirmación
usar exploratory como calibration
recalcular los pesos de W_eq con C1
eliminar una familia por producir resultados desfavorables
cambiar el plan o un CSV sin invalidar sus hashes
```

Cualquier cambio necesario antes de C1 exige una enmienda explícita, un nuevo freeze y nuevos hashes. Después de iniciar C1 no se retoca el protocolo confirmatorio.
