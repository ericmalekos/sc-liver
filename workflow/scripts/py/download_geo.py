"""Stage 00 — download + verify a dataset's raw data and organize it into per-sample
10x triplet directories (results/00_download/<ds>/raw/<GSM>/), then write a DONE sentinel.

GEO RAW tarballs ship per-GSM files with mixed naming; this groups them by GSM accession
and normalizes filenames to barcodes/features/matrix. Heuristic by design — verify the
layout for a new dataset and adjust the regexes if needed. Not exercised by the fixture.
"""
import hashlib
import os
import re
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import get_logger  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
url = str(snakemake.params.url)  # noqa: F821
accession = str(snakemake.params.accession)  # noqa: F821
md5 = snakemake.params.md5  # noqa: F821
outdir = Path(snakemake.params.outdir)  # noqa: F821
sentinel = Path(snakemake.output.sentinel)  # noqa: F821

outdir.mkdir(parents=True, exist_ok=True)
tar_url = url if url.endswith(".tar") else url.rstrip("/") + f"/{accession}_RAW.tar"
tar_path = outdir / f"{accession}_RAW.tar"

if not tar_path.exists():
    log.info(f"Downloading {tar_url}")
    urllib.request.urlretrieve(tar_url, tar_path)

if md5:
    h = hashlib.md5()
    with open(tar_path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    assert h.hexdigest() == md5, f"md5 mismatch for {tar_path}"
    log.info("md5 verified")

with tarfile.open(tar_path) as tf:
    tf.extractall(outdir)

# organize GSM-prefixed flat files into per-sample triplet dirs
# Preserve v2 vs v3 feature filename so scanpy auto-detects the format correctly
# (genes.tsv = 10x v2 [2 cols]; features.tsv = v3 [3 cols]).
NAME = {
    "barcodes": "barcodes.tsv.gz",
    "features": "features.tsv.gz",
    "genes": "genes.tsv.gz",
    "matrix": "matrix.mtx.gz",
}
gsm_re = re.compile(r"(GSM\d+)")
moved = 0
for f in list(outdir.glob("*")):
    if f.is_dir() or f.name.endswith(".tar"):
        continue
    m = gsm_re.search(f.name)
    if not m:
        continue
    gsm = m.group(1)
    kind = next((k for k in NAME if k in f.name.lower()), None)
    if kind is None:
        continue
    dest = outdir / gsm
    dest.mkdir(exist_ok=True)
    target = dest / NAME[kind]
    shutil.move(str(f), str(target))
    moved += 1

n_samples = sum(1 for d in outdir.iterdir() if d.is_dir())
sentinel.write_text(f"accession={accession}\nsamples={n_samples}\nfiles_organized={moved}\n")
log.info(f"Organized {moved} files into {n_samples} sample dirs under {outdir}")
