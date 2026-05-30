# Containers

Two images back the reproducible runtime; per-rule `envs/*.yaml` remain the source of truth.

| Image | Dockerfile | Purpose | Maps to rules |
|---|---|---|---|
| `*-base` | `Dockerfile.base` | CPU runtime: Snakemake + conda; builds per-rule envs at runtime | all CPU rules |
| `*-gpu`  | `Dockerfile.gpu`  | CUDA runtime with `integrate_gpu` + `cellbender` envs pre-baked | scVI/scANVI, CellBender (only when `gpu.enabled`) |

## Build + convert to Singularity (for SLURM)

```bash
bash containers/build.sh base     # -> containers/base.sif
bash containers/build.sh gpu      # -> containers/gpu.sif  (optional)
```

Then run with `--use-singularity --use-conda` (Snakemake resolves each rule's conda
env *inside* the container, so the analysis pins are identical to a bare `--use-conda` run).
The SLURM profile already binds `/private/groups` (`singularity-args`).
