### Protocolo C1 para Bratu3D y proximidad al fold

#### Estado

```text
DRAFT_NO_FREEZE
```

Este documento define cómo construir el eje de régimen de Bratu3D sin usar el éxito o el fracaso de los métodos comparados para localizar el fold.

#### Principio de independencia

El fold discreto `lambda_fold_h` no se define mediante:

```text
fallo de Newton-Krylov
fallo de Chebyshev
número de iteraciones de un método comparado
región donde un método empieza a ganar
```

La localización se basa en una rama de continuación y en un proxy espectral del Jacobiano.

`apps/calibrate_b6.cpp` usa Newton-Krylov únicamente como corrector numérico de continuación. Si el corrector falla, reduce el paso. Si no logra recuperar la rama por encima de `min_step`, termina con error y no declara un fold.

#### Proxy espectral

Para cada punto aceptado de la rama inferior se estima el autovalor de menor magnitud del Jacobiano simétrico mediante iteración inversa con GMRES precondicionado.

Cada estimación debe informar:

```text
lambda
spectral_gap
eigen_residual
continuation_iters
n_F
n_JVP
n_prec
```

Un punto no se acepta como evidencia de calibración si la estimación espectral no converge o si el residual del par propio no satisface la tolerancia fijada.

#### Proxy del fold

La herramienta puede imprimir `FOLD_PROXY_SPECTRAL` a partir de los últimos puntos de la rama usando la relación local esperada de un saddle-node entre `spectral_gap^2` y la distancia al fold.

Este valor es únicamente un proxy de calibración:

```text
NO es lambda_fold_h congelado
NO define por sí solo rho_lambda confirmatorio
NO sustituye pseudo-arclength si la auditoría final lo exige
```

Antes del freeze, el proxy debe verificarse en más de una malla y con refinamiento del paso de continuación. Si la estimación no es estable, se implementa pseudo-arclength antes de C1.

#### Eje confirmatorio

Una vez fijado y versionado `lambda_fold_h` para cada tamaño:

```text
rho_lambda = lambda / lambda_fold_h
```

La malla confirmatoria se mantiene como:

```text
0.20
0.40
0.60
0.80
0.90
0.95
0.975
0.99
```

La calibración usa valores distintos:

```text
0.30
0.50
0.70
0.85
0.925
0.9625
0.9825
```

No se eliminan regiones porque Newton-Krylov o Chebyshev pierdan en ellas.

#### Proxies explicativos

Además de `rho_lambda`, se guardan como variables explicativas:

```text
spectral_gap
residual del par propio
dificultad lineal
número de aplicaciones del precondicionador
distancia inicial
tamaño
```

El `spectral_gap` explica condicionamiento, pero no reemplaza el eje preregistrado `rho_lambda`.

#### Gate antes del freeze

Bratu3D no queda listo para C1 hasta cerrar:

```text
1  continuación reproducible de la rama inferior
2  estimador espectral con residual controlado
3  estabilidad de lambda_fold_h entre refinamientos de paso
4  comprobación en al menos dos tamaños antes del freeze
5  decisión explícita sobre necesidad de pseudo-arclength
6  mallas rho_lambda de calibración y confirmación separadas
7  archivos CSV y commit exacto registrados
```
