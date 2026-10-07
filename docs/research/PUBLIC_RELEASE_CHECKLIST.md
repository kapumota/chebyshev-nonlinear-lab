### Checklist para hacer público el repositorio

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
- reproducción desde clon limpio PASS
- CUDA marcado correctamente como ejecutado o `NOT_EXECUTED_NO_GPU`

#### Documentación

- README coherente con el código real
- licencia presente
- derechos para publicar y relicenciar el código importado bajo la licencia elegida confirmados
- `CITATION.cff` presente
- provenance presente
- protocolo de reproducción presente
- limitaciones conocidas documentadas
- casos pequeños descritos como tests, no como evidencia de gran escala

#### Seguridad y privacidad

- revisar secretos y credenciales
- revisar rutas personales innecesarias
- revisar archivos de configuración local
- revisar datasets y licencias
- revisar archivos grandes

Comandos útiles:

```bash
git status --short
git diff --check
git ls-files | grep -E '(^|/)(build|build-mpi|build-verify)/' || true
git ls-files | grep -E '\.(zip|tar|tar\.gz|o|a|so|exe)$' || true
```

#### Metadata de publicación

Antes de la primera release pública revisar manualmente:

- autores y colaboradores en `CITATION.cff`
- URL pública del repositorio
- versión
- fecha de release
- licencia aplicable al código propio

#### Decisión

El repositorio solo cambia de privado a público cuando todos los puntos aplicables estan resueltos.
