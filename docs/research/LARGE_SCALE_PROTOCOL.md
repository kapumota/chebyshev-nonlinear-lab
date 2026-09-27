### Protocolo de gran escala

#### Objetivo

Separar explicitamente validacion numerica de evaluacion cientifica de escalabilidad.

#### Clasificacion de tamaños

| Nivel | Dimension aproximada | Funcion |
|---|---:|---|
| Toy | N < 100 | unit tests, golden tests, regresion |
| Small | 100 <= N < 10^3 | validacion inicial |
| Medium | 10^3 <= N < 10^4 | validacion del core dinamico y sparse |
| Large | 10^4 <= N < 10^5 | experimentos cientificos principales |
| Very large | N >= 10^5 | memoria, matrix-free, paralelismo y escalabilidad |

La clasificacion se interpreta segun la estructura del problema. Un sistema disperso y uno denso con el mismo N tienen costes radicalmente distintos.

#### Regla de evidencia

Los casos toy no pueden sustentar afirmaciones de:

- escalabilidad
- eficiencia HPC
- comportamiento de memoria a gran escala
- superioridad frente a solvers consolidados

#### Metricas minimas

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
- memoria maxima
- tiempo de comunicacion
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

Los papers de gran escala no deben comparar exclusivamente implementaciones internas. S0 determinara los solvers externos apropiados y sus configuraciones reproducibles.

#### Falsificacion

Cada protocolo experimental debe incluir condiciones donde el metodo candidato pueda perder frente a los baselines. No se deben seleccionar solamente casos favorables.
