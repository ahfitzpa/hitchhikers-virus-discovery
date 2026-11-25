"""Extract unique SRA run accessions from headerless MMseqs2 hit tables.

Invoked via Snakemake `script:` by `rules/01.centroid_anchored.smk` and
`rules/01.centroid_unanchored.smk`. Parses the `<run_id>_...` prefix from the
first column of `snakemake.input.hits` and writes one accession per line to
`snakemake.output[0]`.
"""

import os
import tempfile
from pathlib import Path

import pandas as pd


def _atomic_write_lines(dest: str | Path, lines: list[str]) -> None:
    """Write lines to a sibling temp file and atomically rename into place."""
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{dest_path.name}.", suffix=".tmp", dir=dest_path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            for line in lines:
                fh.write(f"{line}\n")
        os.replace(tmp_name, dest_path)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


hits_path = Path(snakemake.input.hits)
output_path = Path(snakemake.output[0])

try:
    df = pd.read_csv(
        hits_path,
        sep="\t",
        names=["query", "target", "pident", "alnlen", "evalue"],
        header=None,
    )
    ids = [str(run_id) for run_id in df["query"].astype(str).str.split("_").str[0].dropna().unique()]
except pd.errors.EmptyDataError:
    ids = []

_atomic_write_lines(output_path, ids)
