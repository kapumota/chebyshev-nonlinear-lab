### Protocolo C1 para B5, Chandrasekhar H

#### Estado

```text
DRAFT_NO_FREEZE
```

Este documento fija la interpretación de B5 para Paper 1. No contiene resultados confirmatorios ni autoriza reutilizar señales exploratorias previas.

#### Formulación

`ChandrasekharOp` discretiza sobre `[0,1]` la ecuación H en la parametrización

```text
x(s) = 1 + (c/2) x(s) integral_0^1 [s/(s+t)] x(t) dt
```

con `c` en `[0,1]`. Con cuadratura de Gauss-Legendre dinámica:

```text
F_i(x) = 1 + x_i (A x)_i - x_i
A_ij   = (c/2) w_j mu_i / (mu_i + mu_j)
```

Los nodos `mu_i` y pesos `w_i` proceden de Gauss-Legendre sobre `[0,1]`.

La matriz `A` no se materializa. La acción `A v` se evalúa directamente con memoria `O(N)` y trabajo `O(N^2)`.

#### Trazabilidad con la tesis

El caso histórico de la tesis usa `lambda = 1/4`. En la parametrización anterior corresponde a:

```text
c = 0.5
```

El caso `N = 8, c = 0.5` se conserva como prueba de regresión matemática. No es una instancia confirmatoria de C1.

#### Separación entre exploración, calibración y confirmación

Los valores observados antes del freeze, incluidos:

```text
c = 0.5
c = 0.99
c = 0.999
c = 0.9999
c = 1.0
```

son `EXPLORATORY`. No se reutilizan como evidencia confirmatoria ni se seleccionan por haber favorecido a un método.

La malla confirmatoria propuesta en el preregistro v3 se expresa mediante:

```text
q = 0.75, 1.75, 2.75, 3.75, 4.75
delta_c = 10^(-q)
c = 1 - delta_c
```

La calibración usa valores diferentes:

```text
q_cal = 0.50, 1.50, 2.50, 3.50, 4.50
delta_c_cal = 10^(-q_cal)
c_cal = 1 - delta_c_cal
```

`c = 1.0` queda reservado como stress test de frontera y no forma parte de H1 confirmatoria.

#### Métodos y outcomes

B5 debe comparar, como mínimo:

```text
Newton-Krylov afinado
Chebyshev inexacto
Chebyshev inexacto con recycling
TF-Chebyshev-GMRES reproducido
```

El baseline externo se añade cuando la interfaz correspondiente esté cerrada.

Outcomes principales:

```text
tasa de éxito
W_eq
t_total
```

Métricas mecanísticas:

```text
n_F
n_JVP
n_D2
n_lin
n_prec
desglose de la primera y segunda resolución
```

`n_JVP` no sustituye `W_eq`, `t_total` ni la tasa de éxito.

#### Recycling

El recycling simple se evalúa como una ablación sistemática en B5. No se presenta como novedad algorítmica.

No se usa un umbral post hoc del número de iteraciones lineales para decidir si se incluye. Cualquier variante más fuerte, como deflación o memoria entre iteraciones no lineales, requiere una decisión separada antes del freeze y no puede activarse porque los resultados preliminares sean favorables.

#### Gate de B5 antes del freeze

B5 no queda listo para confirmación hasta cumplir simultáneamente:

```text
1  implementación Gauss-Legendre validada
2  regresión N = 8, c = 0.5 contra la tesis
3  mallas de calibración y confirmación separadas
4  costos primitivos para W_eq fijados
5  tuning separado de la confirmación
6  política de tamaños fijada antes de observar C1
7  scripts y CSV versionados
8  fallos conservados en el análisis
```
