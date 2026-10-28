### R1.3, instrumentación reproducible y política de ejecución

#### Estado

```text
IMPLEMENTADO COMO INFRAESTRUCTURA
SIN EVIDENCIA CONFIRMATORIA
```

R1.3 separa instrumentación, generación de estados iniciales y paralelización. No modifica las hipótesis H1-H4 del preregistro C1.

#### Ejecutables

```text
apps/calibrate_b6.cpp
    continuación de la rama inferior de Bratu3D
    proxy espectral independiente del éxito de los métodos comparados

apps/run_campaign.cpp
    ejecutor de calibración y confirmación
    una fila CSV por método, estado inicial, semilla y repetición

apps/check_openmp_consistency.cpp
    verificación numérica entre uno y varios hilos
    no mide aceleración
```

#### Estados iniciales

`run_campaign` no calcula una solución de referencia mediante uno de los métodos comparados.

Los estados admitidos son:

```text
natural
continuation
coarse
perturbed
```

Para `continuation` y `coarse` se suministra un vector versionado mediante `--start-file`.

Para `perturbed` también se exige un vector de referencia externo. Los niveles preregistrados son:

```text
near = 0.05
mid  = 0.20
far  = 0.50
```

La perturbación se normaliza respecto de la norma del vector de referencia. Los tres niveles se conservan aunque uno de los métodos gane o pierda sistemáticamente.

#### Repeticiones

La configuración por defecto del ejecutor es:

```text
1 warm-up excluido
7 repeticiones medidas
12 semillas científicas para perturbed
```

Las repeticiones de timing no se consideran réplicas científicas. Las semillas generan estados iniciales distintos y sí forman parte del diseño de robustez.

Los valores definitivos se congelan en el preregistro antes de C1 confirmatoria.

#### Datos crudos

El CSV conserva, como mínimo:

```text
phase
problem
family
n
N
param
method
start_kind
distance
seed
success
nonlinear_iters
n_F
n_JVP
n_D2
n_lin
n_lin_failures
n_prec
t_setup
t_F
t_JVP
t_D2
t_prec
t_total
max_rss_kb
residual_rel
threads
repetition
commit
```

No se calcula `W_eq` durante la ejecución. Los costos primitivos se fijan en una fase de calibración separada y `W_eq` se deriva posteriormente de los datos crudos.

#### Alcance de t_total

`t_total` incluye la construcción estática del proveedor de precondicionamiento de cada repetición y la ejecución completa del solver.

`t_setup` incluye:

```text
construcción estática del proveedor
construcción dinámica realizada por provider.make(x)
```

La lectura de archivos de estados iniciales y la construcción del objeto problema quedan fuera del intervalo medido.

`max_rss_kb` es el máximo residente del proceso hasta ese punto y no una asignación aislada por método. Se conserva como métrica diagnóstica y debe interpretarse con esa limitación.

#### OpenMP

Los núcleos CPU de vectores, Bratu3D y multigrid usan reparto estático para tamaños suficientemente grandes. Chandrasekhar paraleliza las filas de la aplicación matrix-free de `A`.

La verificación de consistencia compara uno y varios hilos. No constituye benchmark de escalabilidad.

Para resultados de tiempo se registra el número de hilos. La afinidad, gobernador de CPU, versión del compilador y demás entorno se registran fuera del CSV mediante el script de captura de entorno antes del freeze.

#### Gate de instrumentación

Antes de C1 confirmatoria deben cumplirse:

```text
1  commit exacto registrado en cada fila
2  árbol Git limpio verificado por el script de campaña
3  warm-up y repeticiones congelados
4  afinidad y número de hilos congelados
5  costos primitivos para W_eq congelados
6  estados continuation/coarse versionados
7  fallos conservados en CSV
8  análisis generado exclusivamente desde datos crudos
```

No se incorporan al manuscrito resultados producidos antes de cerrar este gate.
