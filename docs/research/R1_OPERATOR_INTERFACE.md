### Interfaz de operadores matrix-free, R1.1

#### Estado

```text
IMPLEMENTADO como primera entrega, sin validación científica
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

`Counted<Op>` envuelve un operador y cuenta `n_F`, `n_JVP` y `n_D2` sin alterar resultados. Los métodos reportan los mismos contadores y un test verifica que coinciden. `KrylovResult` registra además `n_lin_iters` y `n_lin_failures` para distinguir trabajo lineal de fallos de convergencia. Esto implementa la contabilidad por componentes exigida por el protocolo.

#### Componentes incluidos

```text
operators/operator.hpp          interfaz, contadores, jvp_fd y d2_fd por diferencias finitas
operators/baseline_adapter.hpp  adaptador de IFunction<N> del baseline
linalg/gmres.hpp                GMRES reiniciado sin matriz, precondicionamiento por la derecha
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

#### Limitaciones conocidas

```text
no hay reciclaje de subespacios entre las dos resoluciones de Chebyshev
el único precondicionador es Jacobi, débil para Bratu3D
el producto vector y los núcleos son secuenciales, sin OpenMP, MPI ni CUDA
los vectores son std::vector<double>, sin matrices dispersas explícitas
solo existe la familia B1 del protocolo
la segunda derivada por diferencias finitas usa un paso con escala epsilon^(1/4) sin calibrar
```

El reciclaje de subespacios es el paso que pone a prueba la hipótesis de la línea L1. Sin él, la segunda resolución cuesta una resolución completa.

#### Verificación

```text
CTest, 32 de 32 tests correctos
    21 del baseline
    11 nuevos en tests/test_operators.cpp
```

Los nuevos tests cubren los conceptos, JVP y segunda derivada analíticos frente a diferencias finitas, el adaptador del baseline, GMRES, ambos métodos en Bratu3D, la coherencia de contadores, la regresión frente al baseline y dos casos de robustez de GMRES.

Los dos tests de robustez comprueban que GMRES informa una ruptura numérica ante un operador nulo y que rechaza `restart == 0` sin entrar en un bucle.

#### Prueba de humo, no es evidencia

`apps/smoke_bratu.cpp` ejecuta ambos métodos en Bratu3D. No sigue un preregistro, usa una repetición y el precondicionador de Jacobi. Sirve para comprobar que el código corre y para anticipar riesgos. El ejecutable valida los tamaños de malla y `NUMA_LAMBDA`; una entrada inválida termina con código de error 2.

Resultados exploratorios originales de R1.1 con lambda 1.0 y lambda 6.0, desde u = 0:

```text
N        lambda  NK n_JVP  Chebyshev n_JVP  iteraciones no lineales (ambos)
9261     1.0     97       152              5
27000    1.0     139      221              5
97336    1.0     278      402              5
9261     6.0     109      175              5
27000    6.0     153      239              5
```

En esos casos ambos métodos necesitan 5 iteraciones no lineales y Chebyshev usa más productos Jacobiano-vector porque resuelve dos sistemas. Es consistente con el riesgo identificado en S0 para la línea L1. No permite concluir nada sobre el método porque faltan reciclaje, un precondicionador adecuado, afinación preregistrada y problemas más exigentes. Estos resultados no deben citarse como evidencia confirmatoria en manuscritos.
