"""Stage 02 — ambient RNA (Python path: `none` passthrough or `cellbender`).

decontX / SoupX are handled by the R sibling (qc_ambient.R); they work on the FILTERED
matrices GEO ships. CellBender needs raw (unfiltered) droplets and a GPU to be practical;
when those aren't available it logs and passes through so the pipeline still completes.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import get_logger, read_h5ad, write_h5ad  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
method = snakemake.params.method  # noqa: F821
adata = read_h5ad(snakemake.input.h5ad)  # noqa: F821

if method == "cellbender":
    log.warning(
        "cellbender requested via the Python path. CellBender needs the raw "
        "(unfiltered) droplet matrix and a GPU; for the GEO filtered matrices used here "
        "decontX (R path) is recommended. Passing counts through unchanged."
    )
adata.obs["ambient_fraction"] = 0.0
adata.uns["ambient_method"] = method
write_h5ad(adata, snakemake.output.h5ad)  # noqa: F821
log.info(f"ambient method={method}: wrote {adata.n_obs} cells")
