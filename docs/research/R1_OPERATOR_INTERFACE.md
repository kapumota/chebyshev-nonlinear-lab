### Interfaz de operadores matrix-free, R1.1 y R1.2

#### Estado

```text
IMPLEMENTADO como primera y segunda entrega, sin validación científica
```

Esta entrega crea infraestructura. No contiene afirmaciones de novedad ni de superioridad de ningún método. Las preguntas científicas están en `STATE_OF_THE_ART.md` y los experimentos deben seguir `LARGE_SCALE_PROTOCOL.md`.

#### Motivación

El baseline materializa `hessian(x)` como `Matrix<N, N*N>` con N fijo en compilación. La memoria crece como N^3 dobles, unos 8 GB para N = 10^3 y unos 8 TB para N = 10^4. Esa representación no sirve para las escalas del protocolo, por lo que R1 reemplaza el acceso por operadores.

#### Interfaz

Archivo `include/numa/operators/operator.hpp`, espacio de nombres `numa::op`.

```text
size()            dimensión N
eval(x, out)      out = F(x)
jvp(x, v, out)    out = F'(x) v
d2(x, v, out)     out = F''(x)[v, v]          opcional
```

Los conceptos `NonlinearOperator` y `SecondOrderOperator` verifican estos requisitos en compilación. El llamador dimensiona `out`. Los vectores son `std::vector<double>`.

#### Contabilidad

`Counted<Op>` envuelve un operador y cuenta `n_F`, `n_JVP` y `n_D2` sin alterar resultados. Los métodos reportan los mismos contadores y un test verifica que coinciden. `KrylovResult` registra además `n_lin_iters`, `n_lin_failures` y `n_prec` para distinguir trabajo lineal, fallos de convergencia y aplicaciones del precondicionador. Esto implementa la contabilidad por componentes exigida por el protocolo.

#### Componentes incluidos

```text
operators/operator.hpp          interfaz, contadores, jvp_fd y d2_fd por diferencias finitas
operators/baseline_adapter.hpp  adaptador de IFunction<N> del baseline
linalg/gmres.hpp                GMRES reiniciado sin matriz, precondicionamiento por la derecha,
                                guardado de la base de Arnoldi y segunda resolución reciclada
linalg/poisson_mg.hpp           multigrid geométrico de ciclo V para el Laplaciano 3D
methods/krylov_newton.hpp       newton_krylov y chebyshev_inexact
problems/bratu3d.hpp            familia B1, con JVP y F'' analíticos y precondicionador de Jacobi
```

Convención de los métodos: `J z = F(x)` y el paso es `x - z`. En Chebyshev, `x+ = x - z - (1/2) u` con `J u = F''(x)[z, z]`.

#### Decisiones de diseño

Los métodos nuevos conviven con el baseline y no modifican ninguna interfaz existente.

El baseline de tamaño fijo se conserva como suite de regresión. Un test ejecuta `chebyshev_inexact` sobre el caso 2x2 de la tesis con tolerancias lineales estrictas y verifica que coincide con `chebyshev<2>`.

Los términos forzantes siguen la elección 2 de Eisenstat-Walker con las salvaguardas de Kelley.

Si el paso completo de Chebyshev no reduce el residual, se usa el paso de Newton con retroceso. La salvaguarda se contabiliza en `n_F`.

Los hiperparámetros de `KrylovParams` son valores iniciales. Deben fijarse en el preregistro de cada campaña y no ajustarse después de observar los resultados confirmatorios.

#### Robustez de GMRES

La implementación de `gmres` distingue entre falta de convergencia y ruptura numérica mediante `GmresResult::converged` y `GmresResult::breakdown`.

Se rechaza `restart == 0` para evitar un bucle sin progreso. También se comprueba que la rotación de Givens tenga un denominador positivo y finito y que los pivotes de la sustitución hacia atrás sean no nulos y finitos.

Una resolución lineal que alcanza el límite de iteraciones sin converger incrementa `n_lin_failures`, pero puede entregar la mejor aproximación disponible. Una ruptura numérica en la primera resolución detiene la iteración no lineal.

En `chebyshev_inexact`, si la segunda resolución lineal sufre una ruptura numérica, no se utiliza la corrección de orden superior. El método intenta el paso de Newton con retroceso usando la primera dirección calculada.

#### Reciclaje de subespacio, R1.2

Las dos resoluciones lineales de Chebyshev comparten el mismo Jacobiano en cada iteración no lineal. Tras la primera resolución, GMRES puede conservar la base de Arnoldi del último ciclo. `gmres_recycled` proyecta el segundo lado derecho sobre ese subespacio y continúa con GMRES sobre el residual restante.

La implementación reutiliza únicamente el último ciclo de la primera resolución. No implementa GCRO-DR, deflación mediante vectores de Ritz armónicos ni memoria entre iteraciones no lineales. Un test verifica que la proyección reciclada no introduce productos Jacobiano-vector ocultos.

Si el espacio guardado es inconsistente o la proyección reciclada no puede utilizarse de forma segura, el código vuelve a GMRES sin recycling.

#### Precondicionador multigrid

`PoissonMG` implementa un ciclo V geométrico para el Laplaciano discreto 3D y se usa como precondicionador aproximado de Bratu3D. Un test verifica que construye la jerarquía y que un ciclo V reduce el residual; otro comprueba el comportamiento de GMRES con multigrid al refinar los tamaños ensayados.

La interfaz contabiliza cada aplicación del precondicionador mediante `n_prec`. Este contador debe guardarse junto con `n_F`, `n_JVP`, `n_D2` y `n_lin` en toda campaña experimental.

#### Limitaciones conocidas

```text
el recycling usa solo el último ciclo de la primera resolución
no hay deflación ni memoria de subespacios entre iteraciones no lineales
el multigrid aproxima el Jacobiano mediante su parte de Laplaciano
los núcleos introducidos en R1.2 siguen siendo secuenciales
los vectores son std::vector<double>, sin matrices dispersas explícitas
solo existe la familia B1 del protocolo
la segunda derivada por diferencias finitas usa un paso con escala epsilon^(1/4) sin calibrar
los hiperparámetros de KrylovParams no están afinados para evidencia confirmatoria
```

#### Verificación

```text
CTest, 36 de 36 tests correctos
    21 del baseline
    15 nuevos en tests/test_operators.cpp
```

Los tests nuevos cubren los conceptos, JVP y segunda derivada analíticos frente a diferencias finitas, el adaptador del baseline, GMRES, ambos métodos en Bratu3D, la coherencia de contadores, la regresión frente al baseline, la robustez de GMRES, el multigrid y el recycling.

En particular, los tests de R1.2 verifican que un ciclo V reduce el residual, que GMRES precondicionado con multigrid mantiene acotado el número de iteraciones en los tamaños ensayados, que la segunda resolución reciclada converge sin productos Jacobiano-vector ocultos y que Chebyshev con recycling conserva la contabilidad.

#### Prueba de humo, no es evidencia

`apps/smoke_bratu.cpp` ejecuta Newton-Krylov, Chebyshev inexacto y Chebyshev inexacto con recycling en Bratu3D desde `u = 0`. El precondicionador se selecciona mediante `NUMA_PRECOND=mg` o `NUMA_PRECOND=jacobi`. El ejecutable valida `NUMA_LAMBDA`, `NUMA_PRECOND` y los tamaños de malla. Una entrada inválida termina con código de error 2.

Las siguientes observaciones se reprodujeron durante la integración de R1.2 con `n = 16`, `N = 4096` y `lambda = 1.0`. Se omiten los tiempos porque cada caso usa una sola repetición y no constituye evidencia de rendimiento.

```text
precondicionador  método           iter  n_F  n_JVP  n_D2  n_prec  res_rel
multigrid         NK                  4    5      11     0      15  9.9e-11
multigrid         Chebyshev           4    5      20     4      28  1.0e-10
multigrid         Cheb+reciclaje      4    5      18     4      30  1.0e-10
jacobi            NK                  5    6      74     0      79  4.8e-09
jacobi            Chebyshev           5    6     123     5     133  4.8e-09
jacobi            Cheb+reciclaje      5    6      99     5     112  4.8e-09
```

En este smoke test, multigrid reduce fuertemente el trabajo lineal respecto de Jacobi. Con Jacobi, recycling reduce `n_JVP` de Chebyshev de 123 a 99; con multigrid, la reducción es de 20 a 18. Estas observaciones sólo verifican integración y contabilidad. No siguen afinación preregistrada, no incluyen repeticiones y no deben citarse como evidencia confirmatoria en manuscritos.
