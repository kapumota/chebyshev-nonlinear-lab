### Preregistro C1, borrador reforzado v3

#### Estado

```text
DRAFT_NO_FREEZE
```

Este documento es un borrador. No debe usarse como preregistro confirmatorio hasta cerrar `S0_FINAL`, completar los baselines obligatorios, fijar reglas estadísticas y crear un commit y una etiqueta antes de ejecutar C1.

La ruta `experiments/paper3/` se conserva temporalmente por compatibilidad con la serie histórica de parches. Científicamente, C1 corresponde a Paper 1. La migración de rutas se realizará en un commit separado después de integrar la serie y antes de la liberación pública.

#### Identidad

```text
campaña                 C1
línea S0                L1, mapa de régimen y punto de equilibrio
objetivo de manuscrito  Paper 1
ruta histórica          experiments/paper3
commit del código       por fijar antes del freeze
fecha de freeze         por fijar
tag de preregistro      por fijar
```

#### Pregunta principal

En qué régimen, cuantificado por condicionamiento, distancia inicial, dificultad de las resoluciones lineales, tamaño y costo del precondicionador, Chebyshev inexacto matrix-free reduce el costo completo frente a Newton-Krylov afinado, y en qué régimen Newton-Krylov continúa siendo superior.

#### Claims excluidos

```text
Chebyshev de orden tres como novedad
Chebyshev tensor-free como novedad
Chebyshev matrix-free como novedad
recycling de Krylov como novedad aislada
higher-order AD como novedad aislada
superioridad universal de Chebyshev
```

La contribución candidata es la caracterización reproducible de la frontera de régimen y un modelo de punto de equilibrio capaz de explicar y predecir el cambio de método ganador.

#### Hipótesis

```text
H1  principal
    En el estrato de mayor dificultad definido antes de C1 mediante proxies cuantitativos,
    Chebyshev inexacto alcanza una razón pareada de W_eq respecto de Newton-Krylov no
    mayor que 0.90 sin una pérdida de tasa de éxito superior a 5 puntos porcentuales.

H2  principal
    En el estrato de menor dificultad, Newton-Krylov no es materialmente peor que
    Chebyshev inexacto en W_eq y conserva igual o mayor robustez.

H3  secundaria
    El recycling reduce el trabajo de la segunda resolución respecto de Chebyshev sin
    recycling. No se interpreta como novedad algorítmica.

H4  secundaria
    Un modelo de frontera construido sólo con datos de calibración predice el signo de
    la ventaja en instancias confirmatorias no usadas para ajustarlo.
```

Los umbrales anteriores son provisionales hasta el freeze. Sólo pueden modificarse usando evidencia de repetibilidad y ruido de calibración, nunca resultados confirmatorios.

#### Familias de problemas

```text
B1/B6  Bratu 3D sobre la rama inferior, desde régimen moderado hasta proximidad al fold
B3a    Burgers 2D con solución manufacturada, usado como benchmark controlado
B3b    Burgers o problema no manufacturado tomado de literatura o software público
B5     Chandrasekhar H desde régimen moderado hasta proximidad al límite conservativo
Bext   al menos un benchmark o aplicación pública no manufacturada para validez externa
```

B2 y B4 pueden añadirse antes del freeze únicamente si tienen implementación, precondicionador, calibración y pruebas equivalentes. No se añaden después de observar C1 primaria para mejorar el resultado.

#### Bratu3D, lambda y ramas

`lambda` se trata como parámetro de continuación y bifurcación. `lambda = 1.0` es una instancia válida de smoke o referencia, pero no define por sí sola el diseño confirmatorio.

Para cada tamaño se estima `lambda_fold_h` antes del freeze mediante continuación versionada y un proxy espectral independiente del solver comparado. La variable de régimen es:

```text
rho_lambda = lambda / lambda_fold_h
```

Niveles confirmatorios:

```text
rho_lambda = 0.20, 0.40, 0.60, 0.80, 0.90, 0.95, 0.975, 0.99
```

Niveles de calibración, distintos de los confirmatorios:

```text
rho_lambda_cal = 0.30, 0.50, 0.70, 0.85, 0.925, 0.9625, 0.9825
```

La campaña principal sigue la rama inferior por continuación desde `lambda` pequeño. Cada instancia registra:

```text
lambda
lambda_fold_h
rho_lambda
branch_id
```

El fold no se define por convergencia o fallo de `nk`, `cheb`, `rec` ni de ningún baseline. Las ramas superiores quedan fuera del análisis principal y sólo pueden aparecer como análisis suplementario preregistrado.

Los valores absolutos de `lambda` ya observados en smoke tests o pilotos se etiquetan `EXPLORATORY` y no se utilizan para justificar a posteriori los niveles confirmatorios.

#### Chandrasekhar H

Se define:

```text
delta_c = 1 - c
```

Los valores confirmatorios se generan mediante una regla logarítmica fijada antes de C1:

```text
q = 0.75, 1.75, 2.75, 3.75, 4.75
delta_c = 10^(-q)
c = 1 - delta_c
```

La calibración usa una malla desplazada:

```text
q_cal = 0.50, 1.50, 2.50, 3.50, 4.50
delta_c_cal = 10^(-q_cal)
c_cal = 1 - delta_c_cal
```

Los valores `c = 0.5`, `0.99`, `0.9999` y `1.0` observados previamente se consideran exploratorios. No forman parte de H1 confirmatoria. `c = 1.0` puede reportarse como stress test suplementario.

#### Burgers 2D

`B3a` usa solución manufacturada para verificar exactitud y estudiar un Jacobiano no simétrico bajo condiciones controladas. No se presenta como aplicación real.

Los ejes principales de B3a se expresan mediante variables adimensionales o proxies equivalentes fijados antes del freeze, incluyendo un Reynolds de celda o magnitud equivalente y un indicador de dificultad lineal.

`B3b` se selecciona antes del freeze a partir de literatura, software o datos públicos. Sus parámetros principales deben proceder de la fuente externa y no de una búsqueda donde gane uno de los métodos.

#### Validez externa

C1 no se firma como diseño suficiente para Paper 1 hasta seleccionar `Bext` o `B3b` no manufacturado con:

```text
fuente pública o publicación identificable
parámetros reproducibles
criterio de convergencia comparable
escala suficiente
sin selección posterior basada en el método ganador
```

#### Métodos y baselines

```text
nk
    Newton-Krylov interno afinado con Eisenstat-Walker o forcing equivalente y búsqueda
    lineal documentada.

cheb
    Chebyshev inexacto matrix-free sin recycling.

rec
    Chebyshev inexacto matrix-free con recycling.

tf_cheb
    TF-Chebyshev-GMRES reproducido a partir de Eustaquio, Ribeiro y Dumett, 2018.

external
    PETSc SNESNEWTONLS o alternativa externa equivalente en las familias dispersas
    principales.
```

El baseline externo es obligatorio para la evaluación de competitividad práctica y no puede eliminarse después de observar pilotos.

#### Precondicionamiento

Dentro de una instancia, los métodos internos usan el mismo precondicionador salvo que una ablación explícita estudie este factor.

```text
Bratu3D        multigrid geométrico o variante fijada en calibración
Burgers 2D     por definir y verificar antes del freeze
Chandrasekhar  por definir y verificar antes del freeze
```

Se registran construcción, actualización y aplicaciones del precondicionador.

#### Puntos iniciales

Cada método recibe exactamente la misma condición inicial. Se distinguen:

```text
natural
    condición estándar del benchmark.

continuation
    solución o predictor del parámetro anterior cuando corresponda.

coarse
    solución de malla gruesa prolongada cuando corresponda.

perturbed
    estrés controlado alrededor de una referencia.
```

Para `perturbed` se mantienen tres estratos `near`, `mid` y `far`. Ningún estrato se elimina porque un método gane claramente en pilotos.

Propuesta de magnitudes relativas para el borrador:

```text
rho_x = 0.05, 0.20, 0.50
```

Se usan 12 semillas preregistradas por estrato cuando la perturbación sea estocástica. Las semillas de calibración y confirmación son distintas.

#### Tamaños

Cada familia escalable debe tener al menos tres tamaños confirmatorios.

Propuesta compatible con jerarquías multigrid:

```text
Bratu3D
    n = 31, 63, 127

Burgers 2D
    n = 127, 255, 511
```

El tamaño máximo puede ajustarse antes del freeze por memoria o costo medido, pero no después de observar qué método gana. Toda modificación se documenta como enmienda previa a C1.

B5 utiliza una secuencia propia adecuada a su estructura. El tamaño no se fuerza a coincidir con las PDE dispersas.

#### Calibración separada de confirmación

La calibración:

```text
usa parámetros y semillas distintos de C1
usa el mismo conjunto de instancias para todos los métodos
asigna el mismo número máximo de configuraciones a cada método
registra el costo total del tuning
fija forcing terms, restart, precondicionamiento, safeguards y tolerancias
no produce evidencia confirmatoria
```

La selección de hiperparámetros usa la siguiente prioridad:

```text
1  maximizar tasa de éxito
2  minimizar W_eq entre configuraciones exitosas
3  usar t_total como verificación de competitividad práctica
```

No se optimiza únicamente `n_JVP`.

#### Costos de primitivas y W_eq

Antes de C1 se ejecutan microbenchmarks independientes para estimar:

```text
c_F
c_JVP
c_D2
c_prec
t_setup
```

Se define:

```text
W_eq = (n_F*c_F + n_JVP*c_JVP + n_D2*c_D2 + n_prec*c_prec + t_setup) / c_JVP
```

Los coeficientes se estiman con el mismo procedimiento para todos los métodos de una instancia y se congelan antes de C1.

#### Métricas

Resultados principales:

```text
W_eq
t_total
tasa de éxito
```

Resultados mecanísticos:

```text
n_F
n_JVP
n_D2
n_prec
n_lin
iteraciones no lineales
t_setup
t_F
t_JVP
t_D2
t_prec
memoria máxima
residual final
```

`n_JVP` no se utiliza como sustituto del costo completo.

#### Repeticiones y unidad estadística

```text
unidad científica
    instancia definida por familia, parámetros, tamaño, tipo de condición inicial y semilla.

repeticiones temporales
    1 warm-up excluido y 7 repeticiones medidas por instancia y método en el subconjunto
    de timing preregistrado.

resumen temporal
    mediana e IQR.

conteos deterministas
    una observación por instancia y método.
```

Afinidad de CPU, hardware, frecuencia, compilador y entorno se registran. Las repeticiones temporales no se tratan como réplicas científicas independientes.

#### Análisis estadístico

Las comparaciones principales son pareadas sobre la misma instancia.

Se reporta:

```text
razón W_eq candidato / nk
razón t_total candidato / nk
intervalos bootstrap del 95% sobre instancias o semillas científicas
cambio en tasa de éxito en puntos porcentuales
performance profile o análisis equivalente que conserve fallos
```

Para razones de costo se usa log-ratio en el cálculo y se reporta la razón de vuelta en escala original. La configuración del bootstrap y sus semillas se fijan antes del freeze.

#### Fallos

No se eliminan silenciosamente corridas fallidas. Cada fallo conserva:

```text
método
instancia
motivo
costo consumido hasta fallo
residual final
iteraciones completadas
```

El costo condicionado a éxito se reporta junto con la tasa de éxito y no la sustituye.

#### Ablaciones

Obligatorias:

```text
cheb frente a rec
precondicionador principal frente a uno débil en un subconjunto preregistrado
contribución de D2 y segunda resolución a W_eq
estratos de condicionamiento
fuentes de condición inicial
```

Las ablaciones se diseñan antes de C1 y no se crean sólo para explicar resultados inesperados.

#### Criterios de parada

Propuesta:

```text
convergencia primaria
    ||F(x_k)|| <= 1e-8 * ||F(x_0)||

reporte secundario
    nivel 1e-10 cuando sea alcanzable sin modificar la configuración

máximo
    50 iteraciones no lineales salvo excepción preregistrada por benchmark

fallo
    máximo de iteraciones, residual no finito, paso no finito o estancamiento definido
    antes del freeze
```

#### Reproducibilidad

Cada fila de datos crudos registra al menos:

```text
problema
subfamilia
parámetros del problema
lambda
lambda_fold_h
rho_lambda
c
delta_c
branch_id
N
método
tipo de condición inicial
nivel de distancia
semilla
éxito
motivo de fallo
proxy de condicionamiento
iteraciones no lineales
n_F
n_JVP
n_D2
n_prec
n_lin
t_setup
t_total
memoria máxima
residual final
commit
hardware
```

Las columnas no aplicables se dejan vacías de forma explícita. Las tablas y figuras se regeneran desde CSV crudos mediante scripts versionados.

#### Regla sobre datos exploratorios previos

Resultados obtenidos antes del freeze, incluidos smoke tests de `lambda = 1.0` y exploraciones de B5, se conservan como `EXPLORATORY`. Pueden motivar preguntas, detectar errores y diseñar calibración, pero no se reutilizan como prueba confirmatoria de H1-H4.

#### Gate editorial de Paper 1

Antes de declarar C1 suficiente para manuscrito deben existir:

```text
1  gap bibliográfico explícito y S0_FINAL
2  claim principal único y falsable
3  TF-Chebyshev-GMRES reproducido
4  baseline externo disponible
5  varias familias de problemas
6  al menos un caso externo no manufacturado
7  tuning separado de confirmación
8  al menos tres tamaños donde proceda
9  métricas de costo, robustez y memoria
10 ablaciones preregistradas
11 análisis pareado con intervalos
12 resultados negativos conservados
13 scripts, CSV, commit, entorno y hardware reproducibles
```

#### Decisiones abiertas antes del freeze

```text
1  confirmar o ajustar los umbrales de H1-H4 usando sólo ruido de calibración
2  algoritmo exacto para lambda_fold_h y proxy espectral de Bratu3D
3  proxy adimensional y precondicionador de B3
4  precondicionador de B5
5  selección exacta de B3b o Bext
6  reproducción verificable de TF-Chebyshev-GMRES
7  versión y configuración del baseline externo
8  presupuesto final de tuning
9  definición exacta de estancamiento
10 subconjunto exacto de timing si el producto cartesiano completo resulta prohibitivo
```

#### Freeze

C1 sólo pasa de `DRAFT_NO_FREEZE` a preregistro confirmatorio cuando todas las decisiones anteriores están resueltas, el árbol está limpio y se crea un commit y una etiqueta específicos.

#### Enmiendas

Sin enmiendas.
