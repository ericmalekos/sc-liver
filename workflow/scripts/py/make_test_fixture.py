#!/usr/bin/env python
"""Generate the tiny SYNTHETIC fixture used for CI / local end-to-end runs.

Writes per-sample 10x MEX triplets (matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz)
under test/data/<sample_id>/ for the samples in test/samples.test.tsv. The data are
synthetic but structured so the full pipeline produces non-trivial, biologically
sensible results:

  * each cell expresses canonical markers of its assigned liver cell type;
  * fibrotic samples (fibrosis_axis >= 2) expand the disease states
    (activated-HSC/myofibroblast, scar-associated macrophage, scar endothelium) and
    upregulate ECM / fibrosis genes -> real cirrhotic-vs-healthy & F2+ DE signal;
  * ligand/receptor genes (SPP1, PDGFB/PDGFRB, TGFB1, TNFSF12/TNFRSF12A, JAG1/NOTCH3 ...)
    are present so cell-cell communication finds the expected fibrotic axes.

Deterministic (seeded). Run once and commit the output:
    python workflow/scripts/py/make_test_fixture.py
"""
from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp
import yaml

REPO = Path(__file__).resolve().parents[3]
SEED = 0


def load_markers() -> dict:
    with open(REPO / "config/markers_liver.yaml") as fh:
        return yaml.safe_load(fh)


def load_gmt_genes() -> list[str]:
    genes: list[str] = []
    with open(REPO / "config/genesets/fibrosis_core.gmt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes += parts[2:]
    return genes


# extra ligand/receptor genes so LIANA recovers the fibrotic crosstalk
LR_EXTRA = [
    "TGFB1", "TGFBR1", "TGFBR2", "PDGFB", "PDGFRB", "SPP1", "CD44", "ITGAV",
    "TNFSF12", "TNFRSF12A", "JAG1", "DLL4", "NOTCH3", "CCL2", "CCR2", "IL1B", "IL1R1",
]

# per-condition cell-type composition (healthy/low vs fibrotic/high)
COMP_HEALTHY = {
    "Hepatocyte": 0.40, "LSEC": 0.15, "Kupffer_resident": 0.12, "Endothelial_vasc": 0.05,
    "HSC_mesenchyme": 0.05, "Cholangiocyte": 0.05, "T_NK": 0.11, "B_plasma": 0.03,
    "Monocyte_MoMac": 0.04,
}
COMP_FIBROTIC = {
    "Hepatocyte": 0.20, "LSEC": 0.08, "Endothelial_vasc": 0.11, "Kupffer_resident": 0.08,
    "Monocyte_MoMac": 0.15, "HSC_mesenchyme": 0.15, "Portal_fibroblast": 0.05,
    "Cholangiocyte": 0.08, "T_NK": 0.07, "B_plasma": 0.03,
}

# which disease-state signature switches on, by cell type, in fibrotic samples
DISEASE_STATE = {
    "HSC_mesenchyme": "activated_HSC_myofibroblast",
    "Portal_fibroblast": "activated_HSC_myofibroblast",
    "Monocyte_MoMac": "scar_associated_macrophage",
    "Kupffer_resident": "scar_associated_macrophage",
    "Endothelial_vasc": "scar_associated_endothelium",
    "LSEC": "scar_associated_endothelium",
}


def build_gene_universe(markers: dict) -> list[str]:
    genes: list[str] = []
    for grp in ("cell_types", "disease_states", "stromal_discriminators"):
        for vals in markers[grp].values():
            genes += vals
    genes += load_gmt_genes() + LR_EXTRA
    genes = sorted(set(genes))
    # pad with background genes so HVG/normalization behave
    genes += [f"BG{n:04d}" for n in range(200)]
    return genes


def cell_type_markers(markers: dict) -> dict:
    return {ct: list(g) for ct, g in markers["cell_types"].items()}


def simulate_sample(sample_id, axis, modality, rng, gene_index, ct_markers, ds_markers):
    fibrotic = axis >= 2
    comp = COMP_FIBROTIC if fibrotic else COMP_HEALTHY
    n_cells = int(rng.integers(700, 1000))
    n_genes = len(gene_index)

    types = list(comp.keys())
    probs = np.array([comp[t] for t in types], dtype=float)
    probs /= probs.sum()
    assignments = rng.choice(types, size=n_cells, p=probs)

    X = rng.poisson(0.15, size=(n_cells, n_genes)).astype(np.float32)  # ambient-ish baseline
    cell_types = []
    for i, ct in enumerate(assignments):
        cell_types.append(ct)
        # canonical markers high
        for g in ct_markers.get(ct, []):
            if g in gene_index:
                X[i, gene_index[g]] += rng.poisson(12)
        # disease-state program in fibrotic samples (a fraction of cells)
        if fibrotic and ct in DISEASE_STATE and rng.random() < 0.7:
            for g in ds_markers[DISEASE_STATE[ct]]:
                if g in gene_index:
                    X[i, gene_index[g]] += rng.poisson(10)
        # snRNA: slightly sparser, lower mito-like background already low
        if modality == "snRNA":
            X[i] = rng.binomial(X[i].astype(int), 0.8)
    # per-cell library-size jitter
    libfac = rng.uniform(0.7, 1.3, size=(n_cells, 1)).astype(np.float32)
    X = np.rint(X * libfac).astype(np.int64)
    return sp.csr_matrix(X), cell_types


def write_10x(outdir: Path, mat_genes_x_cells: sp.spmatrix, genes: list[str], barcodes: list[str]):
    outdir.mkdir(parents=True, exist_ok=True)
    with gzip.open(outdir / "matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, mat_genes_x_cells, field="integer")
    with gzip.open(outdir / "features.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for g in genes:
            w.writerow([g, g, "Gene Expression"])
    with gzip.open(outdir / "barcodes.tsv.gz", "wt", newline="") as fh:
        for b in barcodes:
            fh.write(b + "\n")


def main() -> int:
    markers = load_markers()
    genes = build_gene_universe(markers)
    gene_index = {g: i for i, g in enumerate(genes)}
    ct_markers = cell_type_markers(markers)
    ds_markers = markers["disease_states"]

    import pandas as pd

    sheet = pd.read_csv(REPO / "test/samples.test.tsv", sep="\t", comment="#")
    for _, row in sheet.iterrows():
        sample_id = row["sample_id"]
        rng = np.random.default_rng(SEED + abs(hash(sample_id)) % 100000)
        mat_cells_x_genes, _cts = simulate_sample(
            sample_id, int(row["fibrosis_axis"]), row["modality"], rng,
            gene_index, ct_markers, ds_markers,
        )
        barcodes = [f"{sample_id}_{i:04d}-1" for i in range(mat_cells_x_genes.shape[0])]
        # 10x convention: features x barcodes
        write_10x(REPO / row["data_path"], mat_cells_x_genes.T.tocsr(), genes, barcodes)
        print(f"  wrote {row['data_path']}: {mat_cells_x_genes.shape[0]} cells x {len(genes)} genes")
    print(f"Fixture written for {len(sheet)} samples; {len(genes)} genes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
