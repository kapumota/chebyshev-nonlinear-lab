### Presupuesto común de tuning para C1

#### Estado

```text
FROZEN_BEFORE_CALIBRATION
```

Este artefacto cierra el número de configuraciones candidatas antes de ejecutar
la calibración competitiva de los métodos internos y PETSc.

#### Presupuesto común

Todos los métodos que reciben tuning competitivo usan exactamente:

```text
4 configuraciones candidatas
start = natural
rtol = 1e-8
max_nl = 50
maxlin = 200
warmup = 1
reps = 7
threads = 1
```

La selección permanece:

```text
success -> W_eq -> t_total
```

`n_JVP` sigue siendo una métrica mecanística y no es objetivo de tuning.

#### Precondicionador por familia

Para evitar que el espacio de búsqueda de precondicionadores favorezca a un
método, el precondicionador queda fijado antes de calibración:

```text
B1    jacobi
B5    none
B3a   mg
Bext  jacobi
```

La misma elección se usa en la slice de tuning interno y PETSc.

#### Cuatro candidatos internos

Los parámetros no explorados quedan fijos:

```text
eta_max      = 0.9
gamma        = 0.9
alpha        = 2.0
eta2_factor  = 10.0
eta2_max     = 0.5
```

El presupuesto 2x2 explora únicamente:

```text
t01  eta_init=0.25  restart=30
t02  eta_init=0.25  restart=50
t03  eta_init=0.50  restart=30
t04  eta_init=0.50  restart=50
```

Newton-Krylov, Chebyshev y Chebyshev con recycling reciben exactamente estos
cuatro candidatos, las mismas instancias y las mismas repeticiones.

#### Cuatro candidatos PETSc

PETSc conserva el 2x2 ya definido:

```text
p01  ksp_rtol=1e-4  restart=30
p02  ksp_rtol=1e-4  restart=50
p03  ksp_rtol=1e-5  restart=30
p04  ksp_rtol=1e-5  restart=50
```

No se pretende que `eta_init` y `ksp_rtol` sean parámetros matemáticamente
idénticos. La igualdad exigida es de presupuesto de búsqueda, dos niveles de
control lineal por dos niveles de restart, cuatro configuraciones en total.

#### Slice de calibración

El tuning competitivo usa solamente `start=natural`. Los starts
`continuation`, `coarse` y `perturbed` pertenecen al mapa científico posterior,
no al espacio usado para elegir hiperparámetros.

#### TF-Chebyshev 2018

`tf_cheb` es un baseline reproducido y congelado mediante su protocolo propio.
No se le añade un espacio de tuning ad hoc para igualar artificialmente el
presupuesto.

#### Prohibiciones

Después de observar calibración no se permite:

```text
aumentar o reducir el número de candidatos
añadir nuevos valores de eta_init
añadir nuevos valores de ksp_rtol
añadir nuevos restart
cambiar el precondicionador por familia
cambiar maxlin, max_nl o rtol para favorecer un método
usar confirmación para seleccionar configuración
```
