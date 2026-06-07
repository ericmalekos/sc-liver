"""Stage 01 — assemble the per-dataset HARMONIZED samplesheet (screening Q1).

Combines per-sample metadata onto a common ordinal fibrosis axis (0-4) plus a binary
F2+/low split, records cell counts, and keeps the original label for provenance.
The axis itself is mapped from heterogeneous source labels by make_samplesheet.py via
config.metadata.fibrosis_mapping; here we validate, derive the binary split, and emit.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, read_h5ad  # noqa: E402

import pandas as pd  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
ds = snakemake.params.ds  # noqa: F821
samples = list(snakemake.params.samples)  # noqa: F821
f2plus = int(snakemake.params.f2plus)  # noqa: F821
h5ads = list(snakemake.input.h5ads)  # noqa: F821
out = snakemake.output.sheet  # noqa: F821

rows = []
for sid, h5 in zip(samples, h5ads):
    a = read_h5ad(h5)
    o = a.obs.iloc[0]
    axis = int(o["fibrosis_axis"])
    rows.append(
        dict(
            sample_id=sid,
            dataset=ds,
            donor_id=str(o.get("donor_id", sid)),
            condition=str(o.get("condition", "NA")),
            fibrosis_stage_raw=str(o.get("fibrosis_stage_raw", "NA")),
            fibrosis_axis=axis,
            fibrosis_bin="high" if axis >= f2plus else "low",
            sort_gate=str(o.get("sort_gate", "NA")),
            modality=str(o.get("modality", "NA")),
            n_cells=int(a.n_obs),
        )
    )

df = pd.DataFrame(rows).sort_values(["fibrosis_axis", "sample_id"])
ensure_parent(out)
df.to_csv(out, sep="\t", index=False)
log.info(
    f"[{ds}] harmonized {len(df)} samples | axis range "
    f"{df.fibrosis_axis.min()}-{df.fibrosis_axis.max()} | "
    f"bins: {df.fibrosis_bin.value_counts().to_dict()}"
)
