### Preregistro final candidato de C1

#### Estado

```text
FINAL_CANDIDATE_NOT_FROZEN
```

Este documento reconcilia C1 para Paper 1 después de integrar B3a, Bext, TF2018, PETSc y el análisis estadístico. No autoriza la campaña confirmatoria. El único estado que autoriza C1 es `READY_FOR_C1`, concedido por el gate de readiness después de ratificación del autor.

#### Pregunta científica

Cuándo compensa pagar una corrección Chebyshev inexacta de segundo orden frente a Newton-Krylov afinado, bajo el mismo problema, condición inicial, tolerancia y presupuesto de calibración, y cómo cambia esa frontera con condicionamiento, dificultad Krylov, distancia inicial, tamaño y costo del precondicionador.

#### Claims excluidos

```text
Chebyshev de orden tres como novedad
Chebyshev tensor-free como novedad
Chebyshev matrix-free como novedad
recycling de Krylov como novedad aislada
higher-order AD como novedad aislada
superioridad universal de Chebyshev
```

La contribución candidata es un mapa de régimen reproducible, falsable y auditado, con contabilidad completa, baselines directos y externos, resultados negativos conservados y un modelo de punto de equilibrio.

#### Hipótesis

```text
H1 principal
    En el estrato de mayor dificultad fijado antes de C1, Chebyshev inexacto
    alcanza una razón geométrica pareada W_eq cheb/nk <= 0.90 y la diferencia
    de éxito cheb - nk no es menor que -5 puntos porcentuales.

H2 principal
    En el estrato de menor dificultad fijado antes de C1, Chebyshev no mejora
    materialmente a Newton-Krylov: la razón geométrica pareada W_eq cheb/nk
    es >= 0.90. La robustez de ambos métodos se reporta sin eliminar fallos.

H3 secundaria
    Recycling reduce trabajo lineal respecto de Chebyshev sin recycling.
    Se reporta la fracción pareada de instancias con n_JVP_rec < n_JVP_cheb,
    junto con W_eq, t_total y tasa de éxito. No se interpreta como novedad.

H4 predictiva
    Un modelo de frontera construido exclusivamente con calibration predice
    el signo de log(W_eq_cheb/W_eq_nk) en al menos 75% de las instancias
    confirmatorias donde ambos métodos convergen.
```

H1 y H2 usan el estimador y bootstrap de `C1_STATISTICAL_ANALYSIS.md`. H3 es mecanística. H4 requiere que el modelo y su SHA-256 estén congelados antes de C1.

#### Familias

```text
B1/B6  Bratu 3D sobre rama inferior, desde régimen moderado hasta cercanía al fold
B3a    Burgers 2D manufacturado, benchmark controlado
B5     Chandrasekhar H
Bext   Brusselator público no manufacturado para validez externa
```

B3a no se presenta como aplicación real. Bext no se elige después de observar qué método gana.

#### Bratu3D

Se usa:

```text
rho_lambda = lambda / lambda_fold_h
```

Niveles confirmatorios:

```text
0.20, 0.40, 0.60, 0.80, 0.90, 0.95, 0.975, 0.99
```

Niveles de calibration:

```text
0.30, 0.50, 0.70, 0.85, 0.925, 0.9625, 0.9825
```

`lambda_fold_h` se determina antes de C1 mediante continuación y un proxy independiente del éxito de los solvers comparados. Ningún nivel se elimina porque Newton o Chebyshev pierdan.

#### Chandrasekhar H

Se define:

```text
delta_c = 1 - c
```

Confirmación:

```text
q = 0.75, 1.75, 2.75, 3.75, 4.75
delta_c = 10^(-q)
c = 1 - delta_c
```

Calibration:

```text
q_cal = 0.50, 1.50, 2.50, 3.50, 4.50
delta_c_cal = 10^(-q_cal)
c_cal = 1 - delta_c_cal
```

Los valores exploratorios previos no se reutilizan como evidencia confirmatoria.

#### B3a y Bext

Las mallas numéricas exactas de B3a y Bext se almacenan en el manifest de instancias congelado. El gate exige:

```text
parámetros de calibration y confirmation explícitos
sin solapamiento entre ambas fases
procedencia y regla de construcción documentadas
al menos tres tamaños donde corresponda
ninguna selección basada en el método ganador
SHA-256 del manifest registrado
```

Los valores de smoke no se convierten automáticamente en puntos de calibration o confirmation.

#### Métodos obligatorios

```text
nk
cheb
rec
tf_cheb
petsc
```

`tf_cheb` reproduce el baseline TF-Chebyshev-GMRES de 2018. `petsc` usa SNESNEWTONLS, búsqueda lineal bt, GMRES y MatShell con Jv analítico.

PETSc se calibra sólo con `phase=calibration` y se congela con la regla:

```text
success -> W_eq -> t_total
```

#### Condiciones iniciales

Se conservan cuando corresponda:

```text
natural
continuation
coarse
perturbed near
perturbed mid
perturbed far
```

Para perturbaciones:

```text
rho_x = 0.05, 0.20, 0.50
```

Las semillas de calibration y confirmation son distintas. La lista exacta y los hashes de los starts se congelan en el manifest de starts antes de C1.

#### Costos

```text
W_eq =
    (n_F*c_F
     + n_JVP*c_JVP
     + n_D2*c_D2
     + n_prec*c_prec
     + t_setup)
    / c_JVP
```

`c_F`, `c_JVP`, `c_D2` y `c_prec` proceden exclusivamente de calibration. Los pesos no se recalculan con confirmación. `t_setup` es el costo real de la corrida.

Resultados principales:

```text
W_eq
t_total
tasa de éxito
```

Resultados mecanísticos incluyen `n_F`, `n_JVP`, `n_D2`, `n_prec`, `n_lin`, memoria, residual y fallos.

#### Tuning

Todos los métodos reciben el presupuesto de tuning fijado antes de C1. Las repeticiones temporales se condensan antes de seleccionar.

La selección interna y PETSc conserva:

```text
1  maximizar tasa de éxito
2  minimizar W_eq entre éxitos
3  minimizar t_total entre éxitos
```

No se usa `n_JVP` como objetivo de tuning.

#### Análisis estadístico

El análisis se define en `docs/research/C1_STATISTICAL_ANALYSIS.md` y `scripts/analyze_c1_confirmatory.py`.

Queda fijado:

```text
log-ratios pareados
estimador exp(mean(log-ratio))
bootstrap percentil 95%
10000 remuestreos
semilla maestra 20261029
diferencia de éxito en puntos porcentuales
performance profiles con fallo = +inf
fallos conservados en el denominador
```

Las repeticiones temporales no son réplicas científicas.

#### Prohibiciones después del freeze

Después de `READY_FOR_C1` no se permite:

```text
cambiar H1-H4
cambiar mallas
cambiar semillas o starts
cambiar tuning
cambiar costos de primitivas
cambiar PETSc o TF2018
eliminar near, mid, far
eliminar una familia por resultados desfavorables
cambiar bootstrap, comparaciones o tratamiento de fallos
repetir C1 para escoger la ejecución favorable
```

Un fallo de infraestructura que invalide la campaña completa se documenta como incidente. No autoriza ajustar hipótesis o configuraciones.

#### Ejecución única

C1 se ejecuta una sola vez desde el commit y tag congelados. El manifest de readiness registra el estado previo:

```text
confirmation_started = false
author_ratified = true
C1_READINESS = PASS
```

La campaña produce datos nuevos. Esos datos no pueden retroalimentar el preregistro.

#### Regla editorial

Un resultado donde Newton-Krylov domina amplias regiones no invalida el estudio. El resultado científico es el mapa de régimen, incluida la ausencia de ventaja de Chebyshev donde corresponda.

El manuscrito debe distinguir:

```text
evidencia confirmatoria
análisis secundarios preregistrados
análisis exploratorios posteriores, si existieran
```
