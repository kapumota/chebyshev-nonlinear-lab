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
