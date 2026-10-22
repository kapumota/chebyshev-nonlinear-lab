### Plantilla de preregistro

Copiar a `experiments/<paper>/PREREGISTRATION.md` y fijar en un commit antes de ejecutar la campaña.

#### Identidad

```text
campaña
línea S0
commit del código
fecha
```

#### Hipótesis

Una hipótesis por línea, formulada de manera falsable. No presentar una hipótesis como resultado.

#### Diseño

```text
problemas y tamaños
puntos iniciales y semillas
baselines con versión y configuración
criterio de parada
hardware previsto
repeticiones por configuración
```

#### Métrica primaria y umbral

```text
métrica primaria
umbral de éxito
condición de falsación
decisión asociada, GO, REFRAME, MERGE o KILL
```

#### Casos de falsación

Listar al menos un caso donde el método candidato podría perder frente a los baselines.

#### Enmiendas

Registrar cada cambio posterior con fecha, motivo y efecto sobre las comparaciones.
