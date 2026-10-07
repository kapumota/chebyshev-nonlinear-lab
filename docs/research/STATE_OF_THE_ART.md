### Estado del arte

#### Estado

```text
PENDING_S0_AUDIT
```

Este archivo no contiene todavía una afirmación de novedad.

#### Objetivo de S0

Realizar una auditoría adversarial de las líneas candidatas antes de implementar el nuevo core de investigación.

La auditoría debe cubrir como mínimo:

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
- aplicaciones grandes con datasets públicos

#### Salida requerida

Para cada línea candidata:

```text
pregunta científica
claim potencial
prior art directo
prior art cercano
software comparable
benchmarks
escala
hardware
métricas
limitaciones
hueco verificable
```

La decisión final debe ser una de:

```text
GO
REFRAME
MERGE
KILL
```

#### Regla

No implementar una contribución principal de R1 hasta que su diferencia frente al estado del arte quede formulada de manera falsificable y verificable.
