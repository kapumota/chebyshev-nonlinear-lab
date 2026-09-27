### Estado del arte

#### Estado

```text
PENDING_S0_AUDIT
```

Este archivo no contiene todavia una afirmacion de novedad.

#### Objetivo de S0

Realizar una auditoria adversarial de las lineas candidatas antes de implementar el nuevo core de investigacion.

La auditoria debe cubrir como minimo:

- higher-order nonlinear solvers
- Chebyshev y Chebyshev-Halley
- tensor methods
- Tensor-GMRES y Tensor-Krylov
- Newton-Krylov e inexact Newton
- matrix-free nonlinear methods
- directional second derivatives
- automatic differentiation de orden superior
- adaptive y safeguarded nonlinear solvers
- distributed nonlinear solvers
- communication-aware y communication-avoiding methods
- GPU nonlinear solvers
- aplicaciones grandes con datasets publicos

#### Salida requerida

Para cada linea candidata:

```text
pregunta cientifica
claim potencial
prior art directo
prior art cercano
software comparable
benchmarks
escala
hardware
metricas
limitaciones
hueco verificable
```

La decision final debe ser una de:

```text
GO
REFRAME
MERGE
KILL
```

#### Regla

No implementar una contribucion principal de R1 hasta que su diferencia frente al estado del arte quede formulada de manera falsificable y verificable.
