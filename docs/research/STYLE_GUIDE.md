### Guía de estilo

#### Código

- firmas de funciones, tipos y nombres de API en inglés
- comentarios en español
- cadenas visibles al usuario en español
- C++20
- evitar dependencias nuevas sin justificación técnica o científica
- preservar interfaces del baseline salvo que exista una razón documentada para cambiarlas

#### Markdown

Los encabezados comienzan en nivel tres:

```markdown
### Título
#### Subtitulo
##### Tercer nivel
```

Evitar encabezados de nivel uno y dos dentro de la documentación del proyecto.

#### Prosa

- usar texto técnico directo
- preferir comas y puntos a puntuación decorativa
- evitar guiones largos tipográficos
- evitar comillas tipográficas
- no presentar una hipótesis como resultado
- distinguir implementación, evidencia y conclusión

#### Git

Los mensajes de commit se escriben en español.

Cada commit debe representar una unidad lógica verificable.

Evitar `git add .` cuando existan build directories, resultados crudos o archivos temporales sin revisar.

#### Investigación

Una funcionalidad nueva debe vincularse con una pregunta científica o una necesidad de infraestructura demostrable.

No se debe reclamar novedad antes de completar el gate de estado del arte correspondiente.
