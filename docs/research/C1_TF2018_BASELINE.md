### Baseline tensor-free Chebyshev 2018

#### Estado

```text
STRUCTURAL_REPRODUCTION_NOT_FROZEN
```

Paper 1 incluye un comparador directo basado en Eustaquio, Ribeiro y Dumett, 2018, `A new class of root-finding methods in R^n: the inexact tensor-free Chebyshev-Halley class`, Computational and Applied Mathematics 37, 6654-6675, DOI 10.1007/s40314-018-0709-0.

No se reclama como contribución propia el método tensor-free, Chebyshev-Halley, GMRES ni las dos resoluciones lineales.

#### Variante implementada

Se implementa la especialización `alpha = 0`, correspondiente al caso Chebyshev tensor-free. La convención local resuelve `J z = F` y usa `d1 = -z`.

La curvatura se aproxima mediante:

```text
C(x)d1 = [J(x + h d1)d1 - J(x)d1] / h
```

con:

```text
h = sqrt(eps_machine) * ||x|| / ||d1||
```

acotado numéricamente al intervalo `[2 eps_machine, 0.5]`.

Las tolerancias relativas de las dos resoluciones siguen la forma publicada:

```text
eta1 = min(1e-8, ||F||^2)
eta2 = min(1e-8, ||F||)
```

con piso `2 eps_machine` para evitar tolerancias inferiores a la resolución aritmética.

El segundo sistema para `alpha = 0` es:

```text
J d2 = -0.5 C(x)d1
```

y la actualización es:

```text
x_next = x + d1 + d2
```

#### Alcance de la reproducción

Esta implementación es una reproducción estructural verificable de las ecuaciones publicadas, no una afirmación de equivalencia bit a bit con el software usado por los autores. Para la comparación de gran escala se permite aplicar el mismo precondicionador declarado para los otros métodos; esa decisión debe registrarse en la campaña.

`n_JVP` incluye los dos productos Jacobiano-vector usados para aproximar `C(x)d1`. `n_D2` permanece en cero porque el baseline no solicita la segunda derivada analítica.

#### Gate antes de C1

1. prueba de curvatura contra una D2 analítica;
2. convergencia en al menos un caso pequeño;
3. configuración exacta de GMRES y precondicionador versionada;
4. inclusión como método fijo `tf2018` en los CSV;
5. no incluir `tf2018` en una búsqueda de hiperparámetros diseñada para los métodos internos salvo protocolo explícito y simétrico.
