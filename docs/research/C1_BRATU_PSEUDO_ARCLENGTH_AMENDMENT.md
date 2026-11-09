### Enmienda técnica C1 para pseudo-arclength de Bratu

#### Estado

```text
PREDEFINED_BEFORE_REAL_PSEUDO_ARCLENGTH_RUN
```

Esta enmienda se activa por un fallo de infraestructura del procedimiento de
continuación natural, no por el resultado de `nk`, `cheb`, `rec`, `tf_cheb`
o `petsc`.

#### Evidencia que activa la contingencia

Commit del intento:

```text
af88ccd6a8c1401da1eab0ab8684d9e7e84c2005
```

El corrector de continuación natural para `n=31`, corrida `coarse`, redujo el
paso hasta quedar por debajo de `min_step=1e-4` y terminó con código 3.

Hashes preservados del intento:

```text
a0b4cede378593df5ced402c99aac98f79fed4186788dc6450676fbfcb463e85  bratu_n31_coarse.csv
afccac5aa30a9391290a57d3b765825e39e701730c0c9f6092226dab577b27fc  bratu_n31_coarse.stdout.txt
90f4f95d23c9c70600fe66fe9558e710b3ab73075e8ccf6b09ce221de529144f  bratu_n31_coarse.stderr.txt
```

El último punto natural aceptado todavía tenía `spectral_gap=0.739055`, por
encima de `spectral_stop=0.05`. Por tanto, el proxy diagnóstico no puede
convertirse en `lambda_fold_h`.

#### Regla del fallback

`calibrate_b6_arclength` repite la continuación natural con exactamente los
mismos parámetros del protocolo. Sólo si la continuación natural vuelve a
agotar `min_step` activa pseudo-arclength.

La inicialización pseudo-arclength usa automáticamente los dos últimos estados
naturales aceptados. No se seleccionan manualmente puntos favorables.

El estado aumentado es:

```text
y = (x, lambda)
```

La tangente se aproxima por secante y se normaliza con:

```text
<y,z>_w =
    dot(x_y, x_z) / N
    + lambda_y * lambda_z
```

El corrector resuelve matrix-free el sistema bordeado:

```text
F(x, lambda) = 0

<t, y - y_predictor>_w = 0
```

La primera implementación técnica intentó aplicar `newton_krylov` directamente
al sistema aumentado. Un smoke con `n=7`, sin valor científico, activó el
fallback pero agotó `arc_min_step`. Antes de ejecutar los tamaños reales se
sustituye ese corrector por una eliminación de Schur escalar que reutiliza
exclusivamente el Jacobiano matrix-free de Bratu, GMRES y `PoissonMG`.

En cada iteración del corrector se resuelven:

```text
J a = -F
J b = exp(x)

delta_lambda = (-g - <t_x,a>/N)
               / (t_lambda + <t_x,b>/N)

delta_x = a + b * delta_lambda
```

donde:

```text
g = <t, y - y_predictor>_w
```

La actualización usa búsqueda lineal sobre la norma combinada:

```text
sqrt((||F||_2 / sqrt(N))^2 + g^2)
```

con tolerancia fija `1e-10`. Los solves lineales conservan
`linear_rtol=1e-10`, reinicio 100 y máximo 2000. Esta corrección modifica el
corrector numérico, no el criterio de fold ni los gates estadísticos.

#### Regla para permanecer en la rama inferior

C1 no define el fold mediante el cruce de rama. Si un corrector pseudo-arclength
produce:

```text
lambda_candidate <= lambda_actual
```

o un autovalor espectral no positivo, el punto se rechaza y el paso de
pseudo-arclength se divide por dos.

Se fija:

```text
arc_min_step = 1e-6
arc_max_attempts = 80
```

Si se alcanza `arc_min_step` o el máximo de intentos sin satisfacer el gate:

```text
NO FREEZE
```

No se relaja el protocolo.

#### Gate espectral y freeze

Se mantienen sin cambios:

```text
spectral_stop <= 0.05
eigen_residual <= 1e-7
R2 de gap^2 >= 0.98
discrepancia coarse/fine <= 0.25%
```

El script `freeze_bratu_folds.py` sigue siendo el único componente que puede
producir `C1_BRATU_FOLDS`.

#### Separación respecto de los métodos comparados

Pseudo-arclength se usa exclusivamente para definir el eje físico-numérico
`rho_lambda`. No interviene en la selección de hiperparámetros ni en la
comparación confirmatoria entre métodos.
