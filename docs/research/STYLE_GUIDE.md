### Guia de estilo

#### Codigo

- firmas de funciones, tipos y nombres de API en ingles
- comentarios en español
- cadenas visibles al usuario en español
- C++20
- evitar dependencias nuevas sin justificacion tecnica o cientifica
- preservar interfaces del baseline salvo que exista una razon documentada para cambiarlas

#### Markdown

Los encabezados comienzan en nivel tres:

```markdown
### Titulo
#### Subtitulo
##### Tercer nivel
```

Evitar encabezados de nivel uno y dos dentro de la documentacion del proyecto.

#### Prosa

- usar texto tecnico directo
- preferir comas y puntos a puntuacion decorativa
- evitar guiones largos tipograficos
- evitar comillas tipograficas
- no presentar una hipotesis como resultado
- distinguir implementacion, evidencia y conclusion

#### Git

Los mensajes de commit se escriben en español.

Cada commit debe representar una unidad logica verificable.

Evitar `git add .` cuando existan build directories, resultados crudos o archivos temporales sin revisar.

#### Investigacion

Una funcionalidad nueva debe vincularse con una pregunta cientifica o una necesidad de infraestructura demostrable.

No se debe reclamar novedad antes de completar el gate de estado del arte correspondiente.
