### Roadmap de investigación

#### Principio

El baseline pequeño se conserva para validación. El desarrollo posterior debe responder a preguntas científicas de gran escala y no a la acumulación de nuevas demos elementales.

#### R0, baseline reproducible

Objetivo:

```text
source-v0.0.0
    -> correcciones verificadas
    -> build reproducible
    -> MPI integrado
    -> documentación coherente
    -> baseline-v0.1.0
```

#### S0, auditoría adversarial del estado del arte

Antes de implementar R1 se auditará literatura, software y benchmarks relacionados.

Cada línea candidata recibe una decisión:

```text
GO
REFRAME
MERGE
KILL
```

No se considera obligatorio producir cuatro papers.

#### R1, large-scale core

Objetivos principales:

- dimensión dinámica
- vectores y matrices dispersas
- interfaz de operadores
- solver lineal iterativo
- precondicionamiento
- medición de memoria y coste

El core existente de tamaño fijo se conserva como baseline y suite de regresión.

#### R2, benchmark escalable

Incorporar al menos una familia de problemas cuyo tamaño pueda controlarse sistemáticamente desde aproximadamente 10^3 hasta 10^5 variables o más, sujeto a sparsity y memoria.

#### R3, aplicación pública de escala significativa

Seleccionar el dominio después de S0. La aplicación debe disponer de datos o casos públicos, baselines externos, escala suficiente y criterios reproducibles.

#### R4, métodos de orden superior a gran escala

Investigar formulaciones que eviten materializar tensores completos y que utilicen cuando proceda:

- productos Jacobiano-vector
- derivadas direccionales de segundo orden
- álgebra lineal inexacta
- estrategias adaptativas
- safeguards verificables

#### Líneas de papers, estado provisional

```text
Paper 1  adaptive large-scale nonlinear methods
Paper 2  distributed large-scale nonlinear methods
Paper 3  GPU and higher-order differentiation
Paper 4  large-scale scientific application
```

Los títulos son marcadores de trabajo. S0 puede fusionar, eliminar o reformular cualquiera de estas líneas.
