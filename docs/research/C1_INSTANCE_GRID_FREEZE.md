### Malla de instancias y plan PETSc de C1

#### Estado

```text
PREDEFINED_BEFORE_CALIBRATION
```

Este artefacto fija antes de ejecutar resultados las mallas que habían quedado abiertas para B3a y Bext y la regla de materialización del plan PETSc. Los smoke previos siguen siendo `exploratory` y no se promueven a calibration ni confirmation.

#### Separación de fases

B1 conserva `rho_lambda_cal = 0.30, 0.50, 0.70, 0.85, 0.925, 0.9625, 0.9825` y confirmation `0.20, 0.40, 0.60, 0.80, 0.90, 0.95, 0.975, 0.99`, con `n = 31, 63, 127`. El valor absoluto se obtiene sólo como `lambda = rho_lambda * lambda_fold_h`, usando folds determinados de forma independiente del solver comparado.

B5 usa `N = 64, 128, 256`. Calibration usa `q = 0.50, 1.50, 2.50, 3.50, 4.50`; confirmation usa `q = 0.75, 1.75, 2.75, 3.75, 4.75`, con `c = 1 - 10^(-q)`.

B3a usa calibration `n = 63, 127, 255`, `nu = 0.08, 0.04, 0.02`, `dt = 0.04, 0.08, 0.16`, `amp = 1`, en factorial 3x3. Confirmation usa `n = 127, 255, 511`, `nu = 0.06, 0.03, 0.015`, `dt = 0.06, 0.12, 0.24`, `amp = 1`, también factorial 3x3. La terna smoke `nu=0.05, dt=0.1, amp=1` no aparece.

Bext conserva `A=1`, `B=3`, `alpha=0.02`. Calibration usa `n = 63, 127, 255` y `dt = 0.004, 0.016, 0.064`. Confirmation usa `n = 127, 255, 511` y `dt = 0.006, 0.024, 0.096`. El smoke `dt=0.01` no aparece.

#### Slice de tuning PETSc

El tuning PETSc usa `start=natural`. Los starts `continuation`, `coarse` y `perturbed` siguen perteneciendo al mapa científico y se congelan por separado. Para mantener equidad, la misma slice de tuning debe usarse cuando se cierre el tuning final de los métodos internos.

El precondicionador queda fijado por familia antes de ejecutar:

```text
B1    jacobi
B5    none
B3a   mg
Bext  jacobi
```

Los cuatro candidatos varían sólo `ksp_rtol` y `restart`:

```text
p01  1e-4  30
p02  1e-4  50
p03  1e-5  30
p04  1e-5  50
```

La selección posterior sigue exclusivamente `success -> W_eq -> t_total`.

#### Materialización

Se requiere un archivo de folds con esquema:

```json
{
  "schema_version": 1,
  "artifact": "C1_BRATU_FOLDS",
  "independent_of_compared_solver": true,
  "folds": {
    "31": {"lambda_fold_h": 0.0},
    "63": {"lambda_fold_h": 0.0},
    "127": {"lambda_fold_h": 0.0}
  }
}
```

Los ceros son sólo marcadores del esquema y no son válidos para ejecutar. Con folds reales:

```bash
python3 scripts/materialize_petsc_calibration_plan.py \
  --bratu-folds experiments/c1/frozen/bratu_folds.json \
  --output experiments/c1/frozen/petsc_calibration_plan.json
```

El `source_commit` por defecto es `HEAD`. El generador rechaza solapamientos entre fases, reutilización de los smoke B3a/Bext, folds faltantes o no positivos y configuraciones PETSc incompatibles.
