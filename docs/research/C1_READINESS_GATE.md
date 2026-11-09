### Gate READY_FOR_C1

#### Propósito

`check_c1_readiness.py` es el último gate antes de ejecutar la campaña confirmatoria.

El gate no ejecuta C1. Sólo verifica que el estado previo es inmutable, trazable y completo.

#### Requisitos

El manifest maestro debe declarar:

```text
schema_version = 1
campaign = C1
status = FROZEN_BEFORE_C1
author_ratified = true
confirmation_started = false
source_commit
preregistration_tag
```

También debe incluir artefactos congelados para:

```text
preregistration
instances
starts
primitive_costs
internal_tuning
petsc_freeze
tf2018
hardware
statistical_protocol
statistical_script
campaign_manifest
```

Cada artefacto contiene una ruta y su SHA-256.

#### Estados internos

El checker exige:

```text
PETSc freeze     FROZEN_CALIBRATION_ONLY
preregistration  FROZEN_BEFORE_C1
```

El manifest PETSc debe corresponder al mismo `source_commit`.

#### Git

Antes de conceder readiness:

```text
source_commit es ancestro de HEAD
árbol Git limpio
tag preregistrado apunta a HEAD
sin cambios de código ejecutable entre source_commit y HEAD
```

`source_commit` identifica la versión del código usada para la calibración.
El commit apuntado por el tag de preregistro puede ser posterior porque incorpora
planes, resultados de calibration, manifests y freezes producidos antes de C1.

No se permite `dirty` ni modificar código ejecutable entre ambos commits. El
checker rechaza cambios en CMake, aplicaciones, bibliotecas, scripts, tests y
otras rutas de código protegidas.

#### Campaña no iniciada

El gate rechaza:

```text
confirmation_started = true
archivos *.csv dentro de experiments/c1/confirmation
```

La ausencia de CSV confirmatorios es una barrera operativa para evitar que el freeze se firme después de mirar resultados.

#### Resultado

Sólo un PASS explícito autoriza el siguiente paso:

```text
C1_READINESS=PASS
READY_FOR_C1
```
