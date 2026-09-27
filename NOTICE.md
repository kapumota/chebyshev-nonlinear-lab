### Avisos y dependencias externas

Este repositorio utiliza o descarga dependencias de terceros durante la configuracion y compilacion.

#### Catch2

Se utiliza para tests unitarios y de regresion mediante CMake `FetchContent`. Catch2 conserva su propia licencia y copyright.

#### Google Benchmark

Se utiliza para microbenchmarks mediante CMake `FetchContent`. Google Benchmark conserva su propia licencia y copyright.

#### MPI, OpenMP y CUDA

Las implementaciones y toolchains correspondientes pertenecen a sus respectivos proyectos y proveedores. Su presencia en el flujo de compilacion no modifica sus licencias.

#### Datos externos

Los futuros datasets descargados en `data/external/` deben conservar su fuente, version, licencia y checksum. No se incorporaran al repositorio archivos externos cuya licencia no permita redistribucion.
