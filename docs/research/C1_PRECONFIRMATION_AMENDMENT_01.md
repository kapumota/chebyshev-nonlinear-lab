### C1, amendment preconfirmatorio 01

Fecha: 9 de noviembre de 2026. Corrección autorizada antes de observar datos confirmatorios.

```
confirmation_started=false
no confirmatory data observed
reason=pre-confirmation operational defect correction
source_commit=70e565b726595697ccc09c7a1d58b3a0e2dc8eff
```

El preregistro original `experiments/c1/frozen/preregistration.json`, commit `26a83529be36d5f702821cccd5f616ab6cff5374`, permanece intacto como registro histórico. También permanecen intactos todos los artefactos originales bajo `experiments/c1/frozen/`. Este documento y los manifests v2 constituyen un registro separado. No se crea ni mueve ningún tag.

#### Alcance exclusivo

Sólo cambian:

1. La cobertura de costos primitivos para los tamaños confirmatorios n=511 ya preregistrados en B3a y Bext. Los costos se obtienen con parámetros de calibration, sin confirmation, interpolación, extrapolación ni préstamo de costos de otras mallas.
2. La serialización canónica de `instance_id`, idéntica entre métodos para una misma instancia del plan. Se aplica al ensamblar el dataset, preservando los CSV raw byte a byte.

No cambian H1-H4, familias, mallas confirmatorias, starts, seeds, métodos, tuning, PETSc, TF2018, W_eq como fórmula, bootstrap, comparaciones, tratamiento de fallos, hardware ni source_commit. No se modifica código científico ni el analizador. Los resultados del suplemento no seleccionan candidatos ni actualizan rankings. Los costos y rankings históricos permanecen intactos.

#### Suplemento de calibration

| Familia | Clave nueva (problem, family, n, N, precond, threads) | Parámetros de calibration | Candidatos | Invocaciones |
|---|---|---|---:|---:|
| B3a | burgers, B3a, 511, 522242, mg, 1 | nu={0.08,0.04,0.02}, dt={0.04,0.08,0.16}, amp=1, factorial 3×3 | 4 | 36 |
| Bext | brusselator, Bext, 511, 1022, jacobi, 1 | dt={0.004,0.016,0.064}, A=1, B=3, alpha=0.02 | 4 | 12 |

Se usan t01–t04 de `C1_TUNING_BUDGET_V1.json`, con verificación contra los comandos del plan de calibration interno original. Todos los comandos llevan phase=calibration, start=natural, methods=nk,cheb,rec, warmup=1, reps=7, threads=1 y el source_commit congelado. Las 48 invocaciones producen 48×3×7=1008 filas medidas. Se agregan 48×3=144 warmups, total 1152 llamadas previstas. No se vuelve a ejecutar ninguna instancia de calibration existente.

`run_c1_cost_supplement.py` conserva por comando stdout, stderr, raw, recibo de argumentos/entorno esperado y ledger. Las unidades completas, incluidos fallos de solver, no se repiten. Ante una interrupción recupera unidades con recibo de salida normal y CSV completo. Los intentos incompletos se archivan como infraestructura y se reanuda su comando entero, sin seleccionar repeticiones o métodos según resultados.

#### Costos v2 y procedencia

Después de ejecutar el suplemento sólo en la workstation congelada, `freeze_c1_costs_v2.py --execute` valida:

- Exactamente las 48 unidades, sus hashes y recibos, 1008 filas, fase calibration y sólo las dos claves n=511 autorizadas.
- Métodos, parámetros, candidatos, configuraciones, siete repeticiones, un hilo, commit y warmup=1 registrado en el comando y recibo. Warmup no es una columna de los CSV existentes, por ello se prueba mediante trazabilidad del ejecutable congelado y sus invocaciones.
- Integridad del CSV original de 6048 filas contra su manifest de ejecución y plan original, sin reejecutarlo.

Combina las 6048 filas originales con las 1008 del suplemento, sin filtrar fallos, y ejecuta **el script original** `scripts/freeze_primitive_costs.py`, sin editarlo ni reimplementar su estimador. Se conserva su mediana de t_primitive/n_primitive, reglas de muestras positivas y mínimo de cinco muestras. El cálculo usa todas las filas admisibles del script, incluidos fallos, sin selección por éxito.

La tabla `experiments/c1/amendment01/primitive_costs_v2.csv` debe tener exactamente 14 claves. Las doce originales deben conservar todas sus celdas idénticas. Las dos nuevas proceden sólo del suplemento. Se guardan CSV agregado suplementario, combinado, hashes de raw/logs/recibos y `primitive_costs_provenance.json`. La readiness reproduce el script y comprueba la igualdad byte a byte de la tabla v2. No basta con que exista un archivo con ese nombre.

La tabla v2 no se entrega ahora, porque sus dos valores nuevos requieren mediciones reales. Tampoco se entrega una tabla vacía o con ceros de sustitución. Los manifests preparados contienen status=PENDING_SUPPLEMENTAL_CALIBRATION y sha256=null para ese archivo inexistente. No son manifests listos para C1. Al congelar datos reales la herramienta los completa con hashes verificables y status=FROZEN_BEFORE_C1_AMENDED. Sólo modifica los manifests v2, nunca los originales.

#### Política de instance_id

La política separada `instance_id_policy.json` fija una identidad desde el plan:

```
C1v2|family|problem|n=integer|parameter_key=exact_binary64_hex|...
```

Las claves de parámetros se ordenan lexicográficamente. Los valores usan `float.hex()`, representación reversible y exacta de binary64, independiente de los formatos decimales de los ejecutores. Incluye todos los parámetros científicos de la instancia del plan, incluso los derivados. No usa el texto raw como fuente de identidad. No cambia los números usados en los comandos ni la columna `param`.

Al convertir raw a canonical sólo están autorizadas:

- method: tf2018 → tf_cheb.
- instance_id: texto del ejecutor → identidad de la misma instancia del plan.

Todas las demás celdas deben permanecer exactamente iguales como cadenas, incluidos formatos numéricos, success, t_total, n_F, n_JVP, n_D2, n_prec, seed, start_kind/start_id, distance, rho_x y repetition. Se verifica equivalencia celda por celda antes de W_eq. El dataset con W_eq es un archivo derivado posterior, calculado por el script existente con costos v2. La canonicalización no añade W_eq ni modifica mediciones.

#### Cadena v2 y autorización de ejecución

Los manifests campaign_manifest_v2, preregistration_v2 y c1_readiness_manifest_v2 enlazan el preregistro original, política, costos y procedencia nuevos. Las referencias cruzadas de campaign y preregistration se actualizan sólo en la finalización local. El plan operacional pasa entonces de pendiente a la versión con hashes efectivos del amendment. No puede actualizarse después del inicio de C1.

`check_c1_readiness_v2.py` requiere costos reales, política congelada, procedencia reproducible, cobertura de toda clave confirmatoria y sólo las 14 claves autorizadas. Las claves originales de tamaños usados únicamente en calibration se mantienen, como exige el registro histórico. No admite tamaños adicionales ajenos a las mallas originales y a los dos suplementos.

También verifica confirmation_started=false, ausencia de archivos confirmatorios, Git local para demostrar que el código protegido sigue idéntico a source_commit, hardware y entorno congelados. No requiere ni realiza ninguna operación remota. Puede aceptar las nuevas herramientas y documentos sin commit, pero nunca cambios staged, unstaged o nuevos archivos en rutas protegidas.

Sólo en la workstation, tras los datos reales y todas las comprobaciones, puede emitir:

```
C1_AMENDED_READINESS=PASS
READY_FOR_C1_AMENDED
```

El preflight de ejecución sigue siendo una barrera posterior para dependencias externas, ejecutables, starts y ausencia de resultados. Readiness v2 no ejecuta C1. La preparación entregada en Work sólo declara:

```
C1_AMENDMENT_PREPARATION=READY
C1_CONFIRMATION=NOT_STARTED
```
