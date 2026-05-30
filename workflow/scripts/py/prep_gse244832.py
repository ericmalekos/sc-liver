#!/usr/bin/env python
"""Adapter for the GSE244832 validation dataset.

GSE244832 ships its snRNA data as ONE combined processed matrix plus a per-cell metadata
table (the disease labels are NOT in GEO's per-sample metadata — only in this bundle). This
splits it into per-sample 10x triplet dirs (so the standard ingest handles it) and emits the
harmonized validation samplesheet rows, mapping the author `condition` onto the fibrosis axis.

    NORMAL -> axis 0 (low) ; NAFL/MASL -> axis 1 (low) ; NASH/MASH -> axis 3 (high, F2+)

Run once after downloading GSE244832_hLIVER_processed_files.tar.gz:
    python workflow/scripts/py/prep_gse244832.py
"""
import csv
import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[3]
SRC = REPO / "results/00_download/gse244832/processed_files"
RAW = REPO / "results/00_download/gse244832/raw"
SHEET = REPO / "config/samples.gse244832.tsv"
AXIS = {"NORMAL": 0, "NAFL": 1, "MASL": 1, "NASH": 3, "MASH": 3}
COND = {"NORMAL": "normal", "NAFL": "MASL", "NASH": "MASH"}


def main():
    genes = pd.read_csv(SRC / "hLIVER_genes.csv", header=None)[0].str.strip('"').values
    meta = pd.read_csv(SRC / "hLIVER_metadata.csv")
    meta.columns = ["barcode"] + list(meta.columns[1:])
    samples = meta["orig.ident"].astype(str)
    conditions = meta["condition"].astype(str).str.strip('"')
    barcodes = meta["barcode"].astype(str).values
    print(f"Loading matrix ({len(genes)} genes x {len(meta)} cells)...")
    # the bundle's *.mtx.gz is actually plain-text MatrixMarket (misnamed) -> detect & open right
    mtx_path = SRC / "hLIVER_counts.mtx.gz"
    with open(mtx_path, "rb") as fh:
        is_gz = fh.read(2) == b"\x1f\x8b"
    stream = gzip.open(mtx_path, "rb") if is_gz else open(mtx_path, "rb")
    M = scipy.io.mmread(stream).tocsc()  # genes x cells
    stream.close()
    assert M.shape == (len(genes), len(meta)), f"shape {M.shape} != ({len(genes)},{len(meta)})"

    rows = []
    for sid in sorted(samples.unique()):
        cols = np.where(samples.values == sid)[0]
        cond_raw = conditions.values[cols[0]]
        axis = AXIS.get(cond_raw.upper(), 0)
        sub = M[:, cols]
        d = RAW / sid
        d.mkdir(parents=True, exist_ok=True)
        with gzip.open(d / "matrix.mtx.gz", "wb") as fh:
            scipy.io.mmwrite(fh, sub, field="integer")
        with gzip.open(d / "genes.tsv.gz", "wt", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            for g in genes:
                w.writerow([g, g])
        with gzip.open(d / "barcodes.tsv.gz", "wt", newline="") as fh:
            for c in cols:
                fh.write(barcodes[c] + "\n")
        rows.append(dict(
            sample_id=sid, dataset="gse244832", condition=COND.get(cond_raw.upper(), cond_raw),
            fibrosis_stage_raw=cond_raw, fibrosis_axis=axis, sort_gate="total",
            modality="snRNA", donor_id=sid, data_path=f"results/00_download/gse244832/raw/{sid}",
        ))
        print(f"  {sid}: {cond_raw}->axis{axis}, {len(cols)} cells")

    df = pd.DataFrame(rows)
    df.to_csv(SHEET, sep="\t", index=False)
    print(f"\nWrote {len(df)} validation samples to {SHEET}")
    print(df.groupby(["condition", "fibrosis_axis"]).size())


if __name__ == "__main__":
    main()
