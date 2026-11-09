### Protocolo de gran escala

#### Objetivo

Separar explícitamente validación numérica de evaluación científica de escalabilidad.

#### Clasificación de tamaños

| Nivel | Dimensión aproximada | Función |
|---|---:|---|
| Toy | N < 100 | unit tests, golden tests, regresión |
| Small | 100 <= N < 10^3 | validación inicial |
| Medium | 10^3 <= N < 10^4 | validación del core dinámico y sparse |
| Large | 10^4 <= N < 10^5 | experimentos científicos principales |
| Very large | N >= 10^5 | memoria, matrix-free, paralelismo y escalabilidad |

La clasificación se interpreta según la estructura del problema. Un sistema disperso y uno denso con el mismo N tienen costes radicalmente distintos.

#### Regla de evidencia

Los casos toy no pueden sustentar afirmaciones de:

- escalabilidad
- eficiencia HPC
- comportamiento de memoria a gran escala
- superioridad frente a solvers consolidados

#### Métricas mínimas

Las campañas futuras deben registrar cuando proceda:

- tiempo end-to-end
- iteraciones no lineales
- norma del residual
- norma del paso
- evaluaciones de F
- evaluaciones del Jacobiano
- productos Jacobiano-vector
- evaluaciones direccionales de segundo orden
- iteraciones del solver lineal
- memoria máxima
- tiempo de comunicación
- speedup
- eficiencia paralela
- success rate

#### Strong scaling

Mantener el problema fijo y aumentar recursos.

Registrar al menos:

```text
T_p
S_p = T_1 / T_p
E_p = S_p / p
```

#### Weak scaling

Aumentar el tamaño del problema con los recursos manteniendo aproximadamente constante el trabajo local.

#### Baselines externos

Los papers de gran escala no deben comparar exclusivamente implementaciones internas. S0 determinará los solvers externos apropiados y sus configuraciones reproducibles.

#### Falsificación

Cada protocolo experimental debe incluir condiciones donde el método candidato pueda perder frente a los baselines. No se deben seleccionar solamente casos favorables.

#### Preregistro

Antes de ejecutar cualquier campaña se crea `experiments/<paper>/PREREGISTRATION.md` a partir de `docs/research/PREREGISTRATION_TEMPLATE.md` y se fija en un commit.

El preregistro declara:

```text
hipótesis
métrica primaria
umbral de éxito
condición de falsación
problemas y tamaños
baselines con versión y configuración
puntos iniciales y semillas
criterio de parada
hardware previsto
```

Después de ejecutar no se modifican los umbrales. Cualquier cambio se registra como enmienda con fecha y motivo, y los resultados se reportan frente al preregistro original y frente a la enmienda.

#### Familias de problemas candidatas

Los tamaños se controlan con un parámetro de malla n y se reportan como N. Los niveles siguen la clasificación de tamaños de este documento.

| Id | Problema | Estructura | Función en el protocolo |
|---|---|---|---|
| B1 | Bratu 3D, ecuación de Liouville-Bratu con stencil de 7 puntos | Dispersa, N = n^3 | Familia principal, segunda derivada diagonal y barata |
| B2 | Brusselator estacionario 2D | Dispersa, acoplamiento entre dos campos | Problema usado en la literatura de Halley con diferenciación automática |
| B3 | Burgers 2D | Dispersa, no simétrica | Problema con Jacobiano no simétrico |
| B4 | Reacción-difusión | Dispersa | Segunda familia con términos no lineales distintos a B1 |
| B5 | Chandrasekhar H, ya presente en el repo | Densa, N hasta 10^4 | Caso de falsación, la factorización reutilizable favorece al baseline directo |
| B6 | B1 con parámetro cercano al punto límite, a calibrar | Dispersa, Jacobiano mal condicionado | Caso de falsación por robustez |

Tamaños de referencia para las familias dispersas:

```text
N = 10^4, 3.2 x 10^4, 10^5, 3.2 x 10^5, 10^6
```

El tamaño máximo por campaña depende de la memoria disponible y se declara en el preregistro.

#### Puntos iniciales

Cada problema se evalúa con tres niveles de distancia a la solución de referencia: cercano, intermedio y lejano. Los puntos iniciales se generan con una semilla fija y al menos 10 perturbaciones por nivel. Se reporta la tasa de éxito por nivel y no solo el promedio.

#### Configuración mínima de baselines

Todo baseline se ejecuta con una configuración afinada y documentada. Un baseline sin afinar invalida la comparación.

```text
Newton-Krylov con GMRES
    término forzante de Eisenstat-Walker
    mismo precondicionador y misma tolerancia que el método candidato
    búsqueda lineal documentada
Newton inexacto con tolerancia lineal fija, como referencia secundaria
Broyden de memoria limitada, como referencia secundaria
Solvers externos, versión y opciones registradas
    NonlinearSolve.jl
    PETSc SNES
    SUNDIALS KINSOL
```

La disponibilidad y versión de los solvers externos se verifica durante S0. Hasta entonces la lista es una propuesta.

#### Criterios de parada y éxito

```text
convergencia
    ||F(x_k)|| <= 1e-8 * ||F(x_0)||
    se reporta también el nivel 1e-10
máximo de iteraciones no lineales
    50
fallo
    máximo de iteraciones, valor no finito o estancamiento del residual
```

Una corrida es exitosa solo si alcanza el criterio de convergencia sin valores no finitos. El tiempo de corridas fallidas no se promedia con el de las exitosas.

#### Contabilidad de costo

El costo se reporta por componentes y no en una sola cifra.

```text
n_F      evaluaciones de F
n_JVP    productos Jacobiano-vector
n_D2     evaluaciones direccionales de segundo orden
n_lin    iteraciones lineales totales
n_prec   aplicaciones del precondicionador
t_setup  tiempo de construcción de precondicionador
t_total  tiempo end-to-end
```

Para comparar métodos con distinto uso de derivadas se agrega un costo equivalente en evaluaciones de F, declarando el factor usado para n_JVP y n_D2 según el método de cálculo (diferencias finitas, números duales o modo Taylor). El factor se mide en cada problema y no se asume.

#### Repetición y estadística

Las corridas de tiempo se repiten al menos 5 veces por configuración. Se reportan mediana y rango intercuartílico. La comparación entre métodos se reporta con intervalos y no solo con medias.

#### Campañas

| Campaña | Línea S0 | Pregunta | Métrica primaria |
|---|---|---|---|
| C1 | L1 | Reduce el Chebyshev inexacto sin matriz con Krylov reciclado los productos matriz-vector totales frente a Newton-Krylov afinado | Productos matriz-vector totales y tiempo |
| C2 | L4 | Qué precisión pierde la segunda derivada direccional y cómo afecta al paso | Error relativo de la derivada, iteraciones y tasa de éxito |
| C3 | L3 | Escala el método en GPU y MPI con núcleos fusionados | Speedup, eficiencia, ancho de banda alcanzado |
| C4 | L2 | Acota una política adaptativa la degradación respecto del mejor método fijo | Arrepentimiento, cociente entre el tiempo de la política y el del mejor método fijo |

La campaña C3 queda bloqueada hasta registrar evidencia de hardware GPU en `evidence/`.

#### Condiciones de falsación por campaña

Los umbrales son propuestas y deben fijarse en el preregistro de cada campaña.

```text
C1  la línea pasa a REFRAME si no reduce los productos matriz-vector totales
    en al menos 10 por ciento frente a Newton-Krylov afinado en 3 de 4
    problemas con N de 10^5 o más
C2  la línea se funde con C1 si las diferencias finitas en precisión doble
    no cambian iteraciones ni tasa de éxito frente a números duales
C3  la línea pasa a REFRAME si la fusión de núcleos no reduce el tiempo
    del término de orden tres frente a núcleos separados
C4  la línea se funde con C1 si el arrepentimiento no mejora al de un
    polialgoritmo por defecto de un solver externo
```

Cada campaña incluye al menos un caso donde se espera que el método candidato pierda, por ejemplo B5 y B6.

#### Captura de entorno

Cada corrida guarda un registro de entorno generado con `scripts/capture_env.sh`. El registro incluye CPU, memoria, compilador, CMake, versión de MPI, estado de la GPU cuando exista, y el commit.

#### Formato de resultados

Los datos crudos se guardan en `artifacts/<paper>/` como CSV con una fila por corrida y las columnas:

```text
problema,N,metodo,punto_inicial,semilla,exito,iter_no_lineales,n_F,n_JVP,n_D2,n_lin,n_prec,t_setup,t_total,memoria_max,residual_final,commit
```

Las tablas y figuras de manuscritos se regeneran desde esos archivos con un script versionado. No se editan a mano.
