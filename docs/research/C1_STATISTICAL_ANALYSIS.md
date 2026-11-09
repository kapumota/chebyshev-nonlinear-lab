### Análisis estadístico confirmatorio de C1

#### Estado

```text
PREDEFINED_BEFORE_C1
```

Este protocolo fija el análisis antes de observar resultados confirmatorios. No contiene cifras de C1 y no puede modificarse para explicar resultados inesperados.

#### Unidad científica y pareamiento

La unidad científica se identifica por:

```text
problem
family
n
N
param
instance_id
start_kind
start_id
distance
rho_x
seed
threads
```

`instance_id` conserva parámetros adicionales de cada familia, por ejemplo `nu`, `dt` y amplitud en B3a. `precond` no forma parte de la clave porque pertenece a la configuración congelada de cada método.

Las repeticiones temporales no son réplicas científicas. Para cada método e instancia se condensan primero mediante la mediana de `W_eq` y `t_total`.

#### Métodos y comparaciones

Métodos obligatorios:

```text
nk
cheb
rec
tf_cheb
petsc
```

Comparaciones predefinidas:

```text
cheb / nk
rec / nk
tf_cheb / nk
petsc / nk
rec / cheb
```

Las cuatro primeras sitúan cada candidato o baseline frente a Newton-Krylov. `rec / cheb` evalúa H3.

#### Log-ratios pareados

Para `W_eq` y `t_total`, sobre instancias donde ambos métodos convergen:

```text
L_i = log(M_candidate_i / M_baseline_i)
Delta = mean(L_i)
R = exp(Delta)
```

`R < 1` favorece al candidato y `R > 1` favorece al baseline.

No se imputan costos artificiales a fallos. El costo condicionado a éxito se interpreta junto con la diferencia de éxito y el performance profile.

#### Bootstrap 95%

Se fija:

```text
tipo             percentil pareado
remuestreos      10000
nivel            95%
semilla maestra  20261029
unidad           instancia científica
```

Para log-ratios se remuestrean pares con éxito común. Para éxito se remuestrean pares de indicadores sobre todo el universo de instancias. Las repeticiones temporales nunca se remuestrean como réplicas independientes.

#### Diferencia de éxito

```text
Delta_success =
  100 * (mean(success_candidate) - mean(success_baseline))
```

Se reporta en puntos porcentuales con intervalo bootstrap del 95%.

#### Performance profile conservando fallos

Se construyen perfiles separados para `W_eq` y `t_total`.

Para una instancia `p` y método `s`, si converge:

```text
r_p,s = M_p,s / min_j M_p,j
```

Si falla:

```text
r_p,s = +inf
```

Si todos los métodos fallan, la instancia permanece en el denominador.

Para todo `tau` finito:

```text
rho_s(tau) =
  numero de instancias con r_p,s <= tau
  --------------------------------------
          numero total de instancias
```

Un fallo nunca cuenta como resuelto para un `tau` finito. La meseta final del perfil conserva explícitamente la robustez.

#### Gates

El script aborta si:

```text
phase != confirmation
falta un método obligatorio
aparece un método no preregistrado
se mezclan commits
el commit contiene dirty
hay cobertura desigual entre métodos
hay repeticiones distintas entre métodos
success cambia entre repeticiones
hay repeticiones duplicadas
W_eq o t_total no son positivos y finitos en un éxito
```

No se eliminan `near`, `mid`, `far`, familias, tamaños o semillas porque un método pierda.

#### Salidas

```text
summary_methods.csv
paired_effects.csv
paired_instances.csv
performance_profile.csv
analysis_manifest.json
```

El manifest registra SHA-256 de entrada y salidas, comparaciones, estimador, bootstrap y tratamiento de fallos.

`W_eq`, `t_total` y tasa de éxito son resultados principales. `n_JVP` permanece como métrica mecanística.
