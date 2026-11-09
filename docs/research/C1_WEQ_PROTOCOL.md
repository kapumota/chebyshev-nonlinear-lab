### Protocolo de costos primitivos y W_eq para C1

#### Estado

```text
DRAFT_NO_FREEZE
```

`W_eq` se deriva de datos crudos. No se calcula durante el solver y sus pesos no se ajustan mirando confirmación.

La definición es:

```text
W_eq =
    (n_F*c_F
     + n_JVP*c_JVP
     + n_D2*c_D2
     + n_prec*c_prec
     + t_setup)
    / c_JVP
```

`c_F`, `c_JVP`, `c_D2` y `c_prec` son tiempos medianos por llamada obtenidos únicamente de corridas `phase=calibration`.

#### Flujo

```text
run_campaign --phase calibration
        |
        v
CSV crudo de calibración
        |
        v
scripts/freeze_primitive_costs.py
        |
        v
tabla de costos versionada y frozen
        |
        v
run_campaign --phase confirmation
        |
        v
scripts/compute_weq.py
        |
        v
CSV derivado con W_eq
```

La tabla de costos se une por:

```text
problem
family
n
N
precond
threads
```

`param` se conserva en el CSV crudo y en el CSV derivado, pero deliberadamente no forma parte de la clave de costos. Las mallas de calibración y confirmación usan parámetros distintos, por lo que exigir igualdad de `param` obligaría a recalibrar sobre puntos confirmatorios. Eso está prohibido.

Esta exclusión es válida porque la tabla representa costo computacional por llamada de una primitiva para una familia, tamaño, precondicionador y número de hilos fijos. Antes del freeze se debe comprobar en calibración que esos costos no presentan una dependencia material con `param`. Si la presentan, el protocolo se detiene y se redefine antes de confirmación.

No se permite usar una tabla de otra familia, dimensión, precondicionador o número de hilos.

#### Reglas

```text
1  sólo phase=calibration alimenta la tabla de costos
2  la tabla queda versionada antes de confirmación
3  confirmación no recalibra c_F, c_JVP, c_D2 ni c_prec
4  t_setup se toma de cada corrida porque es costo real de esa ejecución
5  W_eq no reemplaza t_total ni tasa de éxito
6  n_JVP sigue siendo métrica mecanística
7  fallos se conservan aunque W_eq no sea interpretable como éxito
8  c_F y c_JVP requieren un mínimo preregistrado de muestras
9  c_D2 se exige cuando la calibración contiene Chebyshev o recycling
10 c_prec se exige cuando el precondicionador no es none
11 claves duplicadas en una tabla congelada son error
12 costos no finitos o negativos son error
```

Si un costo requerido no tiene suficientes observaciones, el script falla en lugar de imputarlo.

`c_D2 = 0` sólo es admisible para un grupo de calibración sin métodos que usen D2. `c_prec = 0` sólo es admisible cuando `precond=none` y no hubo aplicaciones de precondicionador.

#### Gate antes del freeze

Antes de versionar la tabla definitiva se debe verificar:

```text
calibration solamente
min_samples satisfecho
métodos esperados presentes
variación de costo por param auditada
mismo build y entorno experimental
CSV crudo conservado
SHA-256 del CSV crudo registrado
```

El freeze de costos es anterior a la ejecución confirmatoria. Ningún resultado confirmatorio puede cambiar los pesos de `W_eq`.
