### Caso externo no manufacturado para C1

#### Estado

```text
BEXT_IMPLEMENTED_NOT_FROZEN
```

El caso externo de Paper 1 se deriva del ejemplo público `ex35.cxx` de PETSc para el Brusselator reacción-difusión unidimensional. Su función es aportar una instancia canónica no manufacturada separada de B3a Burgers, que se conserva como benchmark controlado.

#### Modelo fuente

La fuente pública define:

```text
u_t - alpha u_xx = A + u^2 v - (B+1) u
v_t - alpha v_xx = B u - u^2 v
0 < x < 1
A = 1
B = 3
alpha = 1/50
u(x,0) = 1 + sin(2 pi x)
v(x,0) = 3
u(0,t) = u(1,t) = 1
v(0,t) = v(1,t) = 3
```

C1 convierte un único paso temporal en un sistema algebraico no lineal mediante Euler implícito. Esta discretización propia no pretende reproducir bit a bit el integrador IMEX de PETSc; preserva el modelo, parámetros, estado inicial y condiciones de frontera del ejemplo público.

#### Sistema algebraico

Para cada nodo interior se resuelve:

```text
F_u = (u-u_old)/dt - alpha Lap(u) - A - u^2 v + (B+1)u
F_v = (v-v_old)/dt - alpha Lap(v) - B u + u^2 v
```

La implementación proporciona JVP y segunda derivada direccional analíticos. No existe una solución exacta manufacturada ni se usa uno de los métodos comparados para construir la instancia.

#### Papel en Paper 1

```text
B3a   Burgers manufacturado        control y verificación
Bext  Brusselator público          evidencia no manufacturada
```

Bext no sustituye los baselines externos de solver. PETSc SNES sigue siendo un comparador independiente obligatorio antes del freeze confirmatorio.

#### Reglas

1. El estado natural es el estado previo publicado del modelo.
2. Las perturbaciones requieren un `start-file` explícito, igual que el resto de C1.
3. Los fallos se conservan.
4. `W_eq`, `t_total` y tasa de éxito mantienen el mismo papel que en las otras familias.
5. Ningún resultado exploratorio de Bext modifica H1-H4 ni elimina regiones del diseño.
