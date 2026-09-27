### Roadmap de investigacion

#### Principio

El baseline pequeno se conserva para validacion. El desarrollo posterior debe responder a preguntas cientificas de gran escala y no a la acumulacion de nuevas demos elementales.

#### R0, baseline reproducible

Objetivo:

```text
source-v0.0.0
    -> correcciones verificadas
    -> build reproducible
    -> MPI integrado
    -> documentacion coherente
    -> baseline-v0.1.0
```

#### S0, auditoria adversarial del estado del arte

Antes de implementar R1 se auditara literatura, software y benchmarks relacionados.

Cada linea candidata recibe una decision:

```text
GO
REFRAME
MERGE
KILL
```

No se considera obligatorio producir cuatro papers.

#### R1, large-scale core

Objetivos principales:

- dimension dinamica
- vectores y matrices dispersas
- interfaz de operadores
- solver lineal iterativo
- precondicionamiento
- medicion de memoria y coste

El core existente de tamaño fijo se conserva como baseline y suite de regresion.

#### R2, benchmark escalable

Incorporar al menos una familia de problemas cuyo tamaño pueda controlarse sistematicamente desde aproximadamente 10^3 hasta 10^5 variables o mas, sujeto a sparsity y memoria.

#### R3, aplicacion publica de escala significativa

Seleccionar el dominio despues de S0. La aplicacion debe disponer de datos o casos publicos, baselines externos, escala suficiente y criterios reproducibles.

#### R4, metodos de orden superior a gran escala

Investigar formulaciones que eviten materializar tensores completos y que utilicen cuando proceda:

- productos Jacobiano-vector
- derivadas direccionales de segundo orden
- algebra lineal inexacta
- estrategias adaptativas
- safeguards verificables

#### Lineas de papers, estado provisional

```text
Paper 1  adaptive large-scale nonlinear methods
Paper 2  distributed large-scale nonlinear methods
Paper 3  GPU and higher-order differentiation
Paper 4  large-scale scientific application
```

Los titulos son marcadores de trabajo. S0 puede fusionar, eliminar o reformular cualquiera de estas lineas.
