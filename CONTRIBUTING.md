### Contribuir a Chebyshev Nonlinear Lab

El repositorio se desarrolla como software de investigacion reproducible. Cada cambio debe responder a un objetivo tecnico o cientifico identificable y debe poder auditarse mediante Git.

#### Flujo de trabajo

Trabaje desde `main` actualizado y cree una rama corta para una sola fase o gate.

Ejemplos:

```text
r1/large-scale-core
paper1/p1-a-baseline
paper2/p2-b-matrix-free
```

Cada rama debe seguir la secuencia:

```text
implementacion
    -> tests
    -> evidencia
    -> revision
    -> pull request
    -> merge
```

#### Convenciones de codigo

- firmas, tipos y API publicas en ingles
- comentarios y cadenas visibles al usuario en español
- C++20
- evitar cambios de formato no relacionados con la tarea
- no introducir dependencias nuevas sin justificar su necesidad
- no mezclar refactorizaciones estructurales con nuevos resultados cientificos

Consulte `docs/research/STYLE_GUIDE.md`.

#### Tests obligatorios

Para cambios CPU:

```bash
cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DNUMA_BUILD_BENCH=OFF
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Para cambios MPI:

```bash
./scripts/build_mpi.sh
```

Los cambios que afecten CUDA deben indicar explicitamente si fueron ejecutados en hardware NVIDIA. No se debe registrar `PASS` cuando la ruta no fue ejecutada.

#### Evidencia

Los logs crudos locales no se versionan por defecto. Los resultados cientificos que sustenten un paper deben almacenarse siguiendo la estructura de `experiments/` y producir manifiestos, configuraciones y resultados procesados reproducibles.

#### Commits

Los mensajes de commit se escriben en español, en infinitivo y describiendo una sola unidad logica.

Ejemplos:

```text
Agregar operador sparse para el nucleo dinamico
Integrar GMRES como solver lineal iterativo
Registrar protocolo de escalabilidad del Paper 2
```

#### Pull requests

Cada PR debe indicar:

- objetivo
- archivos principales modificados
- tests ejecutados
- evidencia generada
- limitaciones conocidas
- impacto sobre reproducibilidad
- paper o fase relacionada, cuando corresponda
