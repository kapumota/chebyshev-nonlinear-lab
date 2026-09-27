### Procedencia y trazabilidad

#### Paquete de origen

El baseline se importo desde el archivo original:

```text
chebyshev-solver.zip
```

SHA-256:

```text
09bfee9ef1c3a13b0538226436fc8d8fb6fdf1472103d87176dfcae5e430b705
```

El contenido recibido se congelo sin modificaciones en:

```text
source-v0.0.0
```

#### Regla de trazabilidad

El tag `source-v0.0.0` no se modifica, recrea ni mueve.

Las correcciones detectadas durante R0 se registran en commits posteriores y deben conservar evidencia suficiente para distinguir:

```text
estado recibido
    -> defecto observado
    -> correccion
    -> test que demuestra la correccion
```

#### Frontera del baseline

R0 puede corregir:

- errores demostrables
- integracion de build
- tests
- empaquetado
- documentacion
- reproducibilidad

R0 no debe introducir:

- nuevos algoritmos de investigacion
- nuevas afirmaciones de novedad
- campañas experimentales de papers
- optimizaciones elegidas despues de observar resultados cientificos

#### Dependencias externas

Las dependencias obtenidas mediante CMake no se consideran codigo propio del repositorio. Sus licencias y versiones deben conservarse segun corresponda.

Los futuros datasets externos deben registrar como minimo:

```text
fuente
version o fecha
licencia
checksum
procedimiento de descarga
```
