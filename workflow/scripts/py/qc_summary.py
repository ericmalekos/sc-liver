"""Stage 02 — collate per-sample QC metrics into one dataset-level summary table."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger  # noqa: E402

import pandas as pd  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
rows = []
for jf in snakemake.input.metrics:  # noqa: F821
    with open(jf) as fh:
        m = json.load(fh)
    th = m.pop("thresholds", {})
    m["mito_pct_ceiling"] = th.get("mito_pct_ceiling")
    m["n_mads"] = th.get("n_mads")
    rows.append(m)

df = pd.DataFrame(rows).sort_values("sample")
ensure_parent(snakemake.output.tsv)  # noqa: F821
df.to_csv(snakemake.output.tsv, sep="\t", index=False)  # noqa: F821
log.info(f"QC summary: {len(df)} samples, "
         f"{df['n_cells_after'].sum()} cells retained "
         f"(median {df['pct_kept'].median():.0f}% kept)")
