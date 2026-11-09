### Baselines externos de Paper 1

Este directorio conserva manifests pequeños y versionables. Los CSV crudos grandes se almacenan bajo `raw/`, que permanece ignorado por Git salvo su README.

El archivo `petsc_snes_manifest.example.json` es sólo una plantilla y contiene marcadores `POR_FIJAR`; no satisface el gate.

Antes del freeze:

```bash
python3 scripts/validate_external_baseline_manifest.py \
  experiments/paper1/external_baselines/petsc_snes_manifest.json
```

debe terminar con código 0.
