### Checklist para hacer publico el repositorio

#### Git

- rama `main` contiene solamente cambios revisados
- working tree limpio
- `source-v0.0.0` permanece inmutable
- `baseline-v0.1.0` apunta al freeze reproducible
- no existen build directories versionados
- no existen paquetes ZIP o binarios accidentales versionados

#### Reproducibilidad

- CPU/OpenMP 21/21 PASS
- MPI 4/4 PASS
- reproduccion desde clon limpio PASS
- CUDA marcado correctamente como ejecutado o `NOT_EXECUTED_NO_GPU`

#### Documentacion

- README coherente con el codigo real
- licencia presente
- derechos para publicar y relicenciar el codigo importado bajo la licencia elegida confirmados
- `CITATION.cff` presente
- provenance presente
- protocolo de reproduccion presente
- limitaciones conocidas documentadas
- casos pequenos descritos como tests, no como evidencia de gran escala

#### Seguridad y privacidad

- revisar secretos y credenciales
- revisar rutas personales innecesarias
- revisar archivos de configuracion local
- revisar datasets y licencias
- revisar archivos grandes

Comandos utiles:

```bash
git status --short
git diff --check
git ls-files | grep -E '(^|/)(build|build-mpi|build-verify)/' || true
git ls-files | grep -E '\.(zip|tar|tar\.gz|o|a|so|exe)$' || true
```

#### Metadata de publicacion

Antes de la primera release publica revisar manualmente:

- autores y colaboradores en `CITATION.cff`
- URL publica del repositorio
- version
- fecha de release
- licencia aplicable al codigo propio

#### Decision

El repositorio solo cambia de privado a publico cuando todos los puntos aplicables estan resueltos.
