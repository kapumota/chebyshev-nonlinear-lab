### Protocolo de tuning para C1

#### Regla de selección

El tuning usa exclusivamente `phase=calibration`.

La selección por método es lexicográfica:

```text
1  maximizar tasa de éxito
2  minimizar mediana de W_eq entre instancias exitosas
3  minimizar mediana de t_total entre instancias exitosas
```

`n_JVP` no interviene como objetivo de selección.

#### Exactitud y presupuestos fijos

La tolerancia de convergencia y los límites de trabajo no se tunean. Cada ejecución de `select_tuning.py` exige un único valor de:

```text
rtol
max_nl
maxlin
```

El análisis primario usa la tolerancia preregistrada correspondiente. Una tolerancia secundaria se procesa como otro estrato, no como una configuración competidora.

Los parámetros que sí forman la configuración candidata son:

```text
eta_init
eta_max
gamma
alpha
eta2_factor
eta2_max
restart
```

#### Repeticiones e instancias

Las repeticiones de timing se condensan primero mediante la mediana para cada instancia científica.

Una instancia científica queda identificada por:

```text
problem
family
n
N
param
precond
start_kind
start_id
distance
rho_x
seed
threads
```

Cada configuración de un método debe cubrir exactamente las mismas instancias y los mismos identificadores de repetición. Además, todos los métodos deben cubrir el mismo conjunto de instancias y repeticiones.

Todos los métodos deben recibir el mismo número de configuraciones candidatas. Esto controla el presupuesto de tuning y evita favorecer un método mediante una búsqueda más amplia.

#### Procedencia

El tuning final no mezcla commits. El CSV debe corresponder a un único commit limpio.

La opción `--allow-dirty` existe exclusivamente para smoke tests de infraestructura. Una salida producida con un commit `dirty` no se puede congelar como tuning confirmatorio.

#### Ranking auditable

`select_tuning.py` conserva todas las configuraciones candidatas en la salida. Para cada método registra:

```text
rank
selected
success_rate
median_W_eq
median_t_total
```

La fila `selected=1` es la configuración elegida. Conservar el ranking completo permite auditar la elección y evita ocultar configuraciones perdedoras.

#### Hiperparámetros registrados

`run_campaign` registra explícitamente:

```text
rtol
max_nl
eta_init
eta_max
gamma
alpha
eta2_factor
eta2_max
restart
maxlin
```

Las configuraciones elegidas se versionan antes de ejecutar C1 confirmatoria.

#### Separación de métodos

Newton-Krylov, Chebyshev y Chebyshev con recycling pueden seleccionar configuraciones distintas, pero reciben el mismo presupuesto de configuraciones, las mismas instancias de calibración, las mismas repeticiones y los mismos límites de trabajo.

La selección final no usa ninguna instancia confirmatoria.
