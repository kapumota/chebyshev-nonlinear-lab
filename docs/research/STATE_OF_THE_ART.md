### Estado del arte

#### Estado

```text
S0_RAPID_V3_NOT_FINAL
```

Esta versión sustituye `PENDING_S0_AUDIT` por una auditoría adversarial rápida orientada a decidir claims, baselines y diseño experimental. No es todavía `S0_FINAL` y no autoriza una afirmación de novedad absoluta.

#### Qué falta para S0_FINAL

Antes del freeze confirmatorio de C1 deben completarse como mínimo:

```text
1  repetir búsquedas en al menos una base indexada y realizar backward y forward citation chaining
2  leer el texto completo de los competidores directos y registrar notas verificables
3  verificar versiones, opciones, criterios de parada y licencias de los baselines externos
4  cerrar una matriz claim frente a prior art con GO, REFRAME, MERGE o KILL
5  separar claims algorítmicos, de sistemas y de caracterización empírica
6  registrar referencias nuevas aparecidas entre el preregistro y el envío
```

Textos completos prioritarios:

```text
Eustaquio, Ribeiro y Dumett, 2018
Bader, 2005
Feng y Pulliam, 1997
Steihaug y Suleiman, 2013
Tan, Miao, Edelman y Rackauckas, 2025
Wang y Li, 2026
Parks et al., 2006
Yetkin y Ceylan, 2023
Kuether y Steyer, 2024
Pasquale y Markidis, 2026
Shahab, Susanto y Hatzikirou, 2025
```

#### Alcance de esta pasada

La pasada es una revisión rápida adversarial. Combina búsqueda académica, búsqueda web y verificación puntual de metadatos. No debe presentarse como revisión sistemática publicable ni como PRISMA completo.

La ausencia de un trabajo en esta pasada no demuestra ausencia en la literatura. Los huecos siguientes son hipótesis de novedad que deben sobrevivir `S0_FINAL`.

#### Prior art que condiciona Paper 1

```text
Eustaquio, Ribeiro y Dumett, 2018
    Introducen una clase Chebyshev-Halley inexacta tensor-free para sistemas no lineales.
    Cada iteración utiliza dos resoluciones lineales aproximadas y GMRES, y se compara
    con Newton-GMRES. "Chebyshev inexacto + GMRES + tensor-free" no es novedad.

Bader, 2005
    Desarrolla Tensor-Krylov para sistemas no lineales large-scale y estudia problemas
    difíciles con Jacobianos próximos a singularidad. "Segundo orden + Krylov +
    large-scale" tampoco es novedad por sí solo.

Tan, Miao, Edelman y Rackauckas, 2025
    Muestran solvers de orden superior escalables mediante derivadas direccionales con
    Taylor-mode AD. "Higher-order + directional derivatives + large-scale" colisiona.

Wang y Li, 2026
    Estudian Chebyshev adaptativo sin materializar el tensor completo. "Adaptive
    Chebyshev" y "evitar tensor completo" no son claims suficientes.

Parks et al., 2006
    GCRO-DR establece recycling de subespacios de Krylov para secuencias de sistemas.
    Recycling lineal no es una contribución nueva.

Yetkin y Ceylan, 2023
    Integran recycling en Newton-Krylov para power flow de gran escala. No puede
    reclamarse de forma genérica "Newton-Krylov con recycling".

Kuether y Steyer, 2024
    Reutilizan precondicionador y estado Krylov en simulaciones no lineales large-scale.
    La reutilización dentro de un solver no lineal es prior art directo o cercano.

Pasquale y Markidis, 2026
    Estudian la precisión del JVP matrix-free y su efecto sobre JFNK en CPU y GPU.

NonlinearSolve.jl
    Proporciona selección de algoritmos y rutas Jacobian-free para problemas grandes.
    PETSc SNES y otros paquetes constituyen comparadores prácticos obligatorios.
```

#### Claims que no se deben usar

```text
nuevo método Chebyshev de orden tres
primer Chebyshev tensor-free
primer Chebyshev matrix-free
primer solver higher-order con derivadas direccionales
primer Chebyshev adaptativo
primer Newton-Krylov con recycling
recycling de Krylov como novedad aislada
GPU o MPI por sí solos como novedad
```

#### Paper 1, pregunta y contribución defendible

```text
pregunta científica
    Cuándo compensa pagar una corrección Chebyshev inexacta de segundo orden frente
    a Newton-Krylov afinado, bajo el mismo problema, condición inicial, tolerancia,
    precondicionador y presupuesto de calibración.

claim candidato
    Existe una frontera de régimen cuantificable donde cambia el método ganador. La
    frontera depende del condicionamiento, dificultad Krylov, distancia inicial,
    tamaño y costo del precondicionador.

contribución defendible
    Un mapa de régimen matrix-free reproducible y falsable, con baselines directos y
    externos, contabilidad completa, resultados negativos incluidos y un modelo de
    punto de equilibrio que explique la transición entre métodos.
```

El paper no exige que Chebyshev domine globalmente. Un resultado válido puede mostrar una región limitada donde compensa y regiones amplias donde Newton-Krylov es superior.

#### Baselines mínimos de Paper 1

```text
nk
    Newton-Krylov interno afinado con forcing terms y búsqueda lineal documentados.

cheb
    Chebyshev inexacto matrix-free sin recycling.

rec
    Chebyshev inexacto matrix-free con recycling.

tf_cheb
    Reproducción verificable de TF-Chebyshev-GMRES de Eustaquio et al. 2018.

external
    PETSc SNESNEWTONLS o alternativa externa equivalente en las familias dispersas
    donde pueda usarse el mismo residual y criterios de convergencia comparables.
```

El baseline externo no puede eliminarse porque el candidato gane o pierda en pilotos.

#### Capas de evidencia

Paper 1 debe contener tres capas diferenciadas:

```text
A  benchmarks controlados
   Bratu3D, Chandrasekhar H y Burgers 2D para aislar mecanismos numéricos.

B  mapa de régimen
   condicionamiento, distancia inicial, dificultad Krylov, tamaño y precondicionamiento.

C  validez externa
   al menos un benchmark publicado o aplicación pública no manufacturada, con
   parámetros fijados por la fuente externa y no elegidos después de observar C1.
```

Burgers con solución manufacturada se conserva como verificación y benchmark controlado. No se presenta como evidencia de realismo físico. Antes del freeze se selecciona además una instancia Burgers o una aplicación no manufacturada tomada de literatura o software público.

#### Bratu3D, lambda y diagrama de bifurcación

`lambda` es el parámetro de continuación y bifurcación del benchmark canónico Bratu3D. Un valor concreto, como `lambda = 1.0`, define una instancia válida. No convierte el problema en una categoría artificial ni constituye por sí solo evidencia confirmatoria.

La literatura reciente para el cubo unidad sitúa el primer punto de giro continuo cerca de `lambda = 9.9001`. Ese valor es una referencia de consistencia, no un valor impuesto a la discretización.

Para cada malla se estima antes del freeze el primer punto de giro discreto `lambda_fold_h` mediante un procedimiento de continuación versionado e independiente del método comparado. El eje principal es:

```text
rho_lambda = lambda / lambda_fold_h
```

La campaña confirmatoria recorre una curva continua desde régimen moderado hasta proximidad al fold:

```text
rho_lambda = 0.20, 0.40, 0.60, 0.80, 0.90, 0.95, 0.975, 0.99
```

La calibración utiliza niveles intercalados y distintos de los confirmatorios:

```text
rho_lambda_cal = 0.30, 0.50, 0.70, 0.85, 0.925, 0.9625, 0.9825
```

La rama principal es la rama inferior seguida por continuación desde valores pequeños de `lambda`. Cada fila registra `lambda`, `lambda_fold_h`, `rho_lambda` y `branch_id`. Las ramas superiores se reservan para análisis suplementario preregistrado.

El fold no se define por el éxito o fracaso de Newton-Krylov, Chebyshev o recycling. Si una celda próxima al fold falla, el fallo permanece en la tasa de éxito.

#### Chandrasekhar H, parámetro c

Chandrasekhar H es un benchmark clásico asociado con transferencia radiativa. `c` controla el régimen del operador y su proximidad al límite conservativo. El estudio confirmatorio no debe concentrarse únicamente en valores descubiertos durante exploración donde Chebyshev parezca favorable.

Se usa:

```text
delta_c = 1 - c
```

Los valores confirmatorios se generan con una regla fijada antes de C1 y desplazada respecto de los valores exploratorios ya observados. Propuesta:

```text
q = 0.75, 1.75, 2.75, 3.75, 4.75
delta_c = 10^(-q)
c = 1 - delta_c
```

Los valores `c = 0.5`, `0.99`, `0.9999` y `1.0` ya observados durante exploración no forman parte de la prueba confirmatoria principal. Pueden conservarse como evidencia exploratoria o stress tests. `c = 1.0` se trata como caso límite suplementario.

#### Burgers 2D

El benchmark manufacturado se denomina `B3a` y su función es verificar el operador no simétrico y controlar la dificultad Krylov. Sus parámetros se expresan mediante variables adimensionales cuando la discretización lo permita, incluyendo al menos un Reynolds de celda o proxy equivalente fijado antes del freeze.

`B3b` será una instancia no manufacturada tomada de literatura, software o datos públicos. Debe usar parámetros definidos por la fuente externa y no por la respuesta observada de los métodos.

#### Puntos iniciales

No se utiliza una única familia de perturbaciones gaussianas como única evidencia. Cada método recibe exactamente la misma condición inicial y se distinguen cuatro fuentes:

```text
natural
    condición inicial estándar del benchmark.

continuation
    predictor o solución del parámetro previo cuando corresponda.

coarse
    solución de malla gruesa prolongada cuando corresponda.

perturbed
    perturbaciones controladas para estudiar robustez y distancia inicial.
```

Los estratos `near`, `mid` y `far` de la parte perturbada se conservan incluso cuando uno de los métodos domine claramente. No se elimina un estrato por resultados piloto.

#### Condicionamiento y dificultad lineal

Un mapa de régimen debe medir el régimen y no sólo asignar etiquetas. Antes del freeze se fija un proxy reproducible por familia:

```text
Bratu3D
    rho_lambda, distancia al fold discreto y estimación espectral compatible con la
    estructura del Jacobiano.

Chandrasekhar H
    delta_c y un proxy espectral o de cercanía a singularidad.

Burgers 2D
    variable adimensional del régimen, proxy para operador no simétrico e iteraciones
    GMRES por resolución como variable complementaria.
```

#### Modelo de costo

C1 registra costos de primitivas mediante microbenchmarks separados de la campaña confirmatoria:

```text
c_F       costo de una evaluación de F
c_JVP     costo de un producto Jacobiano-vector
c_D2      costo de una derivada direccional de segundo orden
c_prec    costo de una aplicación del precondicionador
t_setup   costo de construcción o actualización del precondicionador
```

Se define:

```text
W_eq = (n_F*c_F + n_JVP*c_JVP + n_D2*c_D2 + n_prec*c_prec + t_setup) / c_JVP
```

`W_eq`, `t_total` y tasa de éxito son resultados principales. `n_JVP` es una métrica mecanística y no sustituye el costo total.

#### Calibración justa

La calibración y la confirmación usan instancias distintas. Cada método recibe:

```text
el mismo conjunto de instancias de calibración
el mismo número máximo de configuraciones evaluadas
la misma regla de selección
el mismo criterio de éxito
```

La regla de selección prioriza:

```text
1  tasa de éxito
2  W_eq entre configuraciones exitosas
3  t_total como verificación práctica
```

No se seleccionan hiperparámetros únicamente por `n_JVP`.

#### Escala y timing

Cada familia escalable debe incluir al menos tres tamaños confirmatorios, preferiblemente separados por al menos un orden de magnitud entre el menor y el mayor cuando memoria y estabilidad lo permitan.

Los tiempos se obtienen con warm-up, afinidad de CPU documentada y varias repeticiones. Las repeticiones temporales no se consideran réplicas científicas independientes.

#### Ablaciones obligatorias

```text
Chebyshev sin recycling frente a Chebyshev con recycling
precondicionamiento débil frente al precondicionador principal donde sea científicamente útil
impacto de D2 y segunda resolución en el modelo de costo
regímenes de condicionamiento y distancia inicial
```

El recycling es un factor experimental y mecanismo de reducción de costo. No es el claim de novedad.

#### Líneas posteriores

```text
Paper 1  GO CONDICIONAL
         mapa de régimen y punto de equilibrio de Chebyshev inexacto matrix-free.

Paper 2  GO CONDICIONAL
         distributed y communication-aware con reutilización del estado lineal.

Paper 3  DEFER
         GPU-resident y kernel fusion sólo si existe contribución arquitectónica propia.

Paper 4  PENDING_S0
         resiliencia frente a fallos numéricos o silent data corruption, sujeta a una
         auditoría independiente.
```

La selección adaptativa y la precisión de derivadas se integran primero como análisis o extensiones de Paper 1, salvo que aparezca una contribución autónoma.

#### Gate de fortaleza para Paper 1

No se ejecuta C1 confirmatoria hasta cumplir:

```text
1  S0_FINAL cerrado
2  TF-Chebyshev-GMRES reproducido o implementado de forma verificable
3  baseline externo disponible en las familias dispersas principales
4  proxies de régimen y condicionamiento fijados
5  calibración separada y justa completada
6  W_eq, t_total y tasa de éxito implementados
7  al menos tres tamaños en cada familia escalable cuando sea factible
8  puntos iniciales naturales y de robustez definidos
9  un benchmark o aplicación externa no manufacturada seleccionada
10 reglas estadísticas, fallos y criterios de parada congelados
11 scripts, CSV crudos y metadatos suficientes para regenerar tablas y figuras
```

Si 2, 3, 4 o 9 no se cumplen, la campaña sigue siendo útil como estudio interno, pero el claim se reduce y no se presenta como mapa de régimen de referencia.

#### Cadenas para S0_FINAL

```text
1  ("Chebyshev" OR "Chebyshev-Halley" OR "Halley") AND
   ("nonlinear systems" OR "systems of nonlinear equations") AND
   ("Newton-Krylov" OR "Jacobian-free" OR "matrix-free" OR "inexact")

2  ("tensor-Krylov" OR "tensor-GMRES" OR "tensor-free") AND
   ("nonlinear equations" OR "nonlinear systems") AND
   ("large-scale" OR sparse)

3  ("Krylov recycling" OR "GCRO-DR" OR "deflated Krylov") AND
   ("Newton-Krylov" OR "nonlinear systems")

4  ("adaptive Chebyshev" OR "safeguarded Chebyshev" OR "adaptive Halley") AND
   ("nonlinear systems" OR "finite element")

5  ("directional derivative" OR "Taylor-mode" OR "automatic differentiation") AND
   ("higher-order nonlinear solver" OR Halley OR Chebyshev)

6  ("communication-avoiding" OR pipelined OR "s-step") AND
   ("Newton-Krylov" OR "nonlinear solver")
```

#### Decisión provisional

```text
L1  GO CONDICIONAL, claim de caracterización de régimen
L2  MERGE o DEFER
L3  GO CONDICIONAL
L4  MERGE
```

La prioridad es cerrar Paper 1 con una comparación que pueda perder, con baselines fuertes, heterogeneidad experimental y un modelo de costo que explique por qué cambia el método ganador.
