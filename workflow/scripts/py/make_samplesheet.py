#!/usr/bin/env python
"""Generate the REAL-run samplesheet by querying GEO (screening Q1 in practice).

For each dataset in the config, fetch the GEO series metadata, parse each GSM's title
and characteristics, infer condition / sort gate / fibrosis label, map the label onto the
harmonized ordinal axis (config.metadata.fibrosis_mapping), and set data_path to the
extracted triplet directory. Heuristic parsing — review the emitted sheet before running.

    python workflow/scripts/py/make_samplesheet.py --config config/config.yaml --out config/samples.tsv
"""
import argparse
import re

import pandas as pd
import yaml


def infer_axis(text: str, mapping: dict):
    t = text.lower()
    # explicit METAVIR / fibrosis stage
    m = re.search(r"\bf\s?([0-4])\b", t)
    if m:
        return int(m.group(1))
    for label, axis in mapping.items():
        if label.lower().replace("_", " ") in t or label.lower() in t:
            return int(axis)
    if "cirrho" in t:
        return 4
    if any(k in t for k in ["healthy", "normal", "control", "uninvolved"]):
        return 0
    return None


def infer_gate(text: str):
    t = text.lower()
    if "cd45+" in t or "cd45 pos" in t or "leuko" in t or "immune" in t:
        return "CD45pos"
    if "cd45-" in t or "cd45 neg" in t or "non-immune" in t:
        return "CD45neg"
    if "pbmc" in t or "blood" in t:
        return "blood"
    return "total"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--destdir", default="results/00_download/_geo_meta")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    mapping = cfg["metadata"]["fibrosis_mapping"]
    import GEOparse  # imported here so the script is usable without GEOparse for --help

    rows = []
    for ds, dcfg in cfg["datasets"].items():
        acc = dcfg["geo_accession"]
        gse = GEOparse.get_GEO(geo=acc, destdir=args.destdir, silent=True)
        for gsm_name, gsm in gse.gsms.items():
            md = gsm.metadata
            organism = " ".join(md.get("organism_ch1", [])).lower()
            if "homo sapiens" not in organism and "human" not in organism:
                continue  # human only (drops mouse samples in GSE136103)
            title = " ".join(md.get("title", []))
            text = " ".join(
                md.get("title", [])
                + md.get("characteristics_ch1", [])
                + md.get("source_name_ch1", [])
            )
            axis = infer_axis(text, mapping)
            if axis is None:
                continue  # skip samples we cannot stage (review separately)
            cond = (
                "cirrhotic" if axis >= 4 else ("healthy" if axis == 0 else "fibrotic")
            )
            # donor = the PATIENT, not the GSM: collapse sort-gate / replicate suffixes
            # (e.g. "Healthy1_Cd45+", "Healthy1_Cd45-A" -> donor "Healthy1") so pseudobulk
            # DE uses one replicate per patient (avoids pseudoreplication).
            m = re.match(r"\s*([A-Za-z]+[ _]?\d+)", title)
            donor = m.group(1).replace(" ", "").replace("_", "") if m else gsm_name
            rows.append(
                dict(
                    sample_id=gsm_name,
                    dataset=ds,
                    condition=cond,
                    fibrosis_stage_raw=text[:80].replace("\t", " "),
                    fibrosis_axis=axis,
                    sort_gate=infer_gate(text),
                    modality=dcfg["modality"],
                    donor_id=donor,
                    data_path=f"results/00_download/{ds}/raw/{gsm_name}",
                )
            )
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(df)} samples to {args.out}")
    print(df.groupby(["dataset", "fibrosis_axis"]).size())


if __name__ == "__main__":
    main()
