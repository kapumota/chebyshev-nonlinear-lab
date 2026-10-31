### Freeze de `lambda_fold_h` para C1

#### Estado

```text
PREDEFINED_BEFORE_FOLD_RUN
```

Este documento fija la regla que convierte la evidencia de continuación y el proxy espectral en `lambda_fold_h`. La regla se congela antes de observar los valores numéricos para `n=31,63,127`.

#### Evidencia requerida

Para cada tamaño se ejecutan dos barridos independientes del mismo ejecutable:

```text
coarse  step = 0.25
fine    step = 0.125
```

Parámetros comunes:

```text
lambda_start   = 0.0
lambda_stop    = 10.2
min_step       = 1e-4
spectral_stop  = 5e-2
eigen_rtol     = 1e-7
linear_rtol    = 1e-10
threads        = 1
```

`lambda_stop` es solamente un límite superior de búsqueda. No se interpreta como fold. Un fallo del corrector no declara el fold.

#### Ajuste local

Para cada corrida se toman los cuatro últimos puntos aceptados y se ajusta:

```text
spectral_gap^2 = a + b * lambda
```

Se exige:

```text
b < 0
R2 >= 0.98
lambda_proxy = -a / b > ultimo lambda aceptado
ultimo spectral_gap <= 0.05
eigen_residual <= 1e-7 en todos los puntos usados
```

El ajuste se vuelve a calcular en Python desde los CSV. No se copia ciegamente el texto `FOLD_PROXY_SPECTRAL` impreso por el ejecutable.

#### Gate de refinamiento

Para cada tamaño:

```text
relative_difference =
    abs(lambda_proxy_fine - lambda_proxy_coarse)
    / abs(lambda_proxy_fine)
```

Se exige:

```text
relative_difference <= 0.0025
```

equivalente a `0.25%`.

Si el gate pasa:

```text
lambda_fold_h = lambda_proxy_fine
```

Si falla cualquiera de los criterios anteriores:

```text
NO FREEZE
implementar pseudo-arclength antes de C1
```

No se permite relajar el umbral después de observar los datos.

#### Ejecución

Después de commitear este protocolo:

```bash
python3 scripts/run_bratu_fold_calibration.py

python3 scripts/freeze_bratu_folds.py   experiments/c1/fold_calibration/bratu_fold_run_manifest.json   experiments/c1/frozen/bratu_folds.json
```

Antes de la corrida real puede ejecutarse:

```bash
python3 scripts/freeze_bratu_folds.py --self-test
```

El self-test sólo valida la lógica del gate con datos sintéticos. No produce evidencia científica.

#### Resultado

El JSON congelado usa:

```text
artifact = C1_BRATU_FOLDS
status = FROZEN_CALIBRATION_ONLY
independent_of_compared_solver = true
```

y contiene, para `31`, `63` y `127`:

```text
lambda_fold_h
coarse_proxy
fine_proxy
relative_difference
fine_r2
hashes de los CSV
source_commit
```

Ese archivo es la única entrada válida para materializar los valores absolutos de `lambda` del plan PETSc.
