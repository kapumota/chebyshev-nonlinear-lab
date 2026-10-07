### Procedencia y trazabilidad

#### Paquete de origen

El baseline se importó desde el archivo original:

```text
chebyshev-solver.zip
```

SHA-256:

```text
09bfee9ef1c3a13b0538226436fc8d8fb6fdf1472103d87176dfcae5e430b705
```

El contenido recibido se congeló sin modificaciones en:

```text
source-v0.0.0
```

#### Regla de trazabilidad

El tag `source-v0.0.0` no se modifica, recrea ni mueve.

Las correcciones detectadas durante R0 se registran en commits posteriores y deben conservar evidencia suficiente para distinguir:

```text
estado recibido
    -> defecto observado
    -> corrección
    -> test que demuestra la corrección
```

#### Frontera del baseline

R0 puede corregir:

- errores demostrables
- integración de build
- tests
- empaquetado
- documentación
- reproducibilidad

R0 no debe introducir:

- nuevos algoritmos de investigación
- nuevas afirmaciones de novedad
- campañas experimentales de papers
- optimizaciones elegidas después de observar resultados científicos

#### Dependencias externas

Las dependencias obtenidas mediante CMake no se consideran código propio del repositorio. Sus licencias y versiones deben conservarse según corresponda.

Los futuros datasets externos deben registrar como mínimo:

```text
fuente
versión o fecha
licencia
checksum
procedimiento de descarga
```
