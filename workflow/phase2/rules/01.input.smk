# =====================================================================
#  inputs.smk: manifest-driven read acquisition
#
#  Handles two sources in one code path:
#     source=sra    -> AWS ODP (or prefetch fallback) + fasterq-dump
#     source=local  -> user's own fastq(.gz) already on disk
#
#  Handles two layouts:
#     layout=paired -> {sample}_1.fastq.gz + {sample}_2.fastq.gz
#     layout=single -> {sample}.fastq.gz            (ONT / PacBio / SE)
#
#  Everything is normalised into ONE directory so downstream rules never
#  need to know where the data came from:
#
#     data/01.reads/{family}/{sample}/{sample}_1.fastq.gz
#     data/01.reads/{family}/{sample}/{sample}_2.fastq.gz     (paired)
#     data/01.reads/{family}/{sample}/{sample}.fastq.gz       (single)
#
#  Downstream rules should call get_reads(wildcards) instead of
#  hard-coding paths.
# =====================================================================

import os
import re
import glob

import pandas as pd
from snakemake.exceptions import WorkflowError


# ---------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------
MANIFEST_PATH = config.get("manifest", "configs/manifest.csv")

READS_DIR = config.get("reads_dir", "data/01.reads")        # normalised reads
SRA_CACHE = config.get("sra_cache", "data/00.shared_sra")   # shared SRA cache

# "symlink" (default, no extra disk) or "copy" (safer on clusters where the
# original data sits on a filesystem the compute nodes cannot see)
LOCAL_LINK_MODE = config.get("local_link_mode", "symlink")

def _find_stage_script():
    """
    workflow.basedir points at the *main* Snakefile, which may sit at the repo
    root or inside workflow/, so probe the usual spots instead of guessing.
    Override with --config stage_script=/path/to/stage_reads.sh
    """
    if config.get("stage_script"):
        return config["stage_script"]
    base = getattr(workflow, "basedir", ".")
    candidates = [
        os.path.join(base, "scripts", "stage_reads.sh"),
        os.path.join(base, "workflow", "scripts", "stage_reads.sh"),
        os.path.join(base, "..", "scripts", "stage_reads.sh"),
        os.path.join("workflow", "scripts", "stage_reads.sh"),
        os.path.join("scripts", "stage_reads.sh"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return os.path.abspath(c)
    raise WorkflowError(
        "[inputs] cannot locate scripts/stage_reads.sh; tried:\n  "
        + "\n  ".join(candidates)
        + "\nSet it explicitly with --config stage_script=/path/to/stage_reads.sh"
    )


STAGE_SCRIPT = _find_stage_script()

# Delete the shared SRA fastq once every family that needs it has been staged.
# Set keep_sra_cache: true to keep them (matches the old behaviour).
KEEP_SRA_CACHE = bool(config.get("keep_sra_cache", False))
_maybe_temp = (lambda p: p) if KEEP_SRA_CACHE else temp

VALID_SOURCES = {"sra", "local"}
VALID_LAYOUTS = {"paired", "single"}
VALID_PLATFORMS = {"illumina", "ont", "pacbio", "bgi", "other"}

_PLATFORM_ALIASES = {
    "nanopore": "ont", "oxford_nanopore": "ont", "oxford nanopore": "ont",
    "minion": "ont", "promethion": "ont", "gridion": "ont",
    "pb": "pacbio", "hifi": "pacbio", "revio": "pacbio", "sequel": "pacbio",
    "novaseq": "illumina", "miseq": "illumina", "nextseq": "illumina",
    "hiseq": "illumina", "mgi": "bgi", "dnbseq": "bgi",
}
_LAYOUT_ALIASES = {
    "pe": "paired", "paired-end": "paired", "paired_end": "paired",
    "se": "single", "single-end": "single", "single_end": "single",
    "long": "single", "nanopore": "single",
}
_SOURCE_ALIASES = {
    "ncbi": "sra", "ena": "sra", "download": "sra", "public": "sra",
    "file": "local", "disk": "local", "own": "local", "in-house": "local",
}

# Columns the loader guarantees exist (missing ones are created empty).
_MANIFEST_COLS = [
    # --- required ---
    "sample_id",       # unique, filesystem-safe; becomes the {sample} wildcard
    "family",          # grouping key used by the rest of the workflow
    "source",          # sra | local   (inferred when blank)
    # --- conditional ---
    "run_id",          # SRR/ERR/DRR accession        (required if source=sra)
    "fq1",             # R1 / ONT reads              (required if source=local)
    "fq2",             # R2                          (required if local+paired)
    # --- strongly recommended ---
    "platform",        # illumina | ont | pacbio     (default illumina)
    "layout",          # paired | single             (inferred when blank)
    "barcode",         # ONT/Illumina barcode or index (provenance only)
    "library_id",      # groups re-sequenced libraries of the same sample
    "batch",           # sequencing run / flowcell ID, for batch-effect checks
    # --- optional metadata, fine to leave blank ---
    "condition",       # treatment / group / phenotype
    "replicate",       # biological or technical replicate number
    "host_species",    # organism / strain
    "collection_date",
    "notes",
    # --- optional QC / tuning hints ---
    "md5_fq1",
    "md5_fq2",
    "adapter",             # override adapter sequence for trimming
    "genome_size",         # e.g. 5m (used by assemblers such as Flye)
    "expected_coverage",
    # --- legacy host/virus-family tracking (from the old
    # download_manifest_parsed.csv schema; get_hostile_ref() in
    # Snakefile_phase2 keys off these). Populated by convert_manifest.py
    # for migrated rows, blank otherwise.
    "bioproject",
    "scientificname",
    "ref_genome_path",
    "scenario",
    "host_status",
]

_SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


# ---------------------------------------------------------------------
# Manifest loading + validation
# ---------------------------------------------------------------------
def _split_paths(spec):
    """
    Expand one manifest path cell into a concrete, sorted list of files.

    Accepts, in any combination separated by ';' or ',':
      * a single file           reads/A_R1.fastq.gz
      * a glob                  reads/A_L00*_R1.fastq.gz
      * a directory             ont/barcode07/          (all fastq inside)
      * several of the above    lane1_R1.fq.gz;lane2_R1.fq.gz

    Multiple files are concatenated at staging time in sorted order, which
    is the normal case for split lanes and for ONT barcode folders.
    """
    out = []
    for part in re.split(r"[;,]", str(spec)):
        part = os.path.expanduser(part.strip())
        if not part:
            continue
        if os.path.isdir(part):
            hits = []
            for pat in ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"):
                hits += glob.glob(os.path.join(part, pat))
            if not hits:
                raise WorkflowError(
                    f"[manifest] directory '{part}' contains no fastq files"
                )
            out += sorted(hits)
        elif any(ch in part for ch in "*?["):
            hits = sorted(glob.glob(part))
            if not hits:
                raise WorkflowError(f"[manifest] glob '{part}' matched nothing")
            out += hits
        else:
            out.append(part)
    return out


def load_manifest(path):
    if not os.path.exists(path):
        raise WorkflowError(
            f"[manifest] not found: {path}\n"
            f"            point at one with --config manifest=/path/to/manifest.csv"
        )

    sep = "\t" if path.endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(path, sep=sep, dtype=str, comment="#").fillna("")
    df.columns = [c.strip().lower() for c in df.columns]

    for col in _MANIFEST_COLS:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].astype(str).str.strip()

    df = df[df["sample_id"] != ""].reset_index(drop=True)
    if df.empty:
        raise WorkflowError(f"[manifest] no usable rows in {path}")

    dupes = df.loc[df["sample_id"].duplicated(), "sample_id"].unique().tolist()
    if dupes:
        raise WorkflowError(f"[manifest] duplicate sample_id: {dupes}")

    problems = []
    for i, row in df.iterrows():
        sid = row["sample_id"]
        ln = i + 2  # header is line 1

        if not _SAFE_ID.match(sid):
            problems.append(
                f"line {ln}: sample_id '{sid}' must be alphanumeric plus . _ -"
            )
        if not row["family"]:
            problems.append(f"line {ln} [{sid}]: 'family' is required")

        # --- source ---
        src = row["source"].lower()
        src = _SOURCE_ALIASES.get(src, src)
        if not src:
            src = "sra" if row["run_id"] else "local"
        if src not in VALID_SOURCES:
            problems.append(
                f"line {ln} [{sid}]: source '{src}' not in {sorted(VALID_SOURCES)}"
            )
        df.at[i, "source"] = src

        # --- platform ---
        plat = row["platform"].lower() or "illumina"
        plat = _PLATFORM_ALIASES.get(plat, plat)
        if plat not in VALID_PLATFORMS:
            problems.append(
                f"line {ln} [{sid}]: platform '{plat}' not in {sorted(VALID_PLATFORMS)}"
            )
        df.at[i, "platform"] = plat

        # --- layout (inferred when blank) ---
        lay = row["layout"].lower()
        lay = _LAYOUT_ALIASES.get(lay, lay)
        if not lay:
            if plat in ("ont", "pacbio"):
                lay = "single"
            elif src == "local":
                lay = "paired" if row["fq2"] else "single"
            else:
                lay = "paired"
        if lay not in VALID_LAYOUTS:
            problems.append(
                f"line {ln} [{sid}]: layout '{lay}' not in {sorted(VALID_LAYOUTS)}"
            )
        df.at[i, "layout"] = lay

        # --- source-specific requirements ---
        if src == "sra":
            if not row["run_id"]:
                problems.append(f"line {ln} [{sid}]: source=sra needs a 'run_id'")
        else:
            if not row["fq1"]:
                problems.append(f"line {ln} [{sid}]: source=local needs 'fq1'")
                continue
            if lay == "paired" and not row["fq2"]:
                problems.append(
                    f"line {ln} [{sid}]: layout=paired needs 'fq2' "
                    f"(use layout=single for ONT / single-end)"
                )
            if lay == "single" and row["fq2"]:
                problems.append(
                    f"line {ln} [{sid}]: layout=single must leave 'fq2' blank"
                )

            # fail fast on unreadable inputs rather than 40 rules later
            for col in ("fq1", "fq2"):
                if not row[col]:
                    continue
                try:
                    files = _split_paths(row[col])
                except WorkflowError as e:
                    problems.append(f"line {ln} [{sid}] {col}: {e}")
                    continue
                for f in files:
                    if not os.path.exists(f):
                        problems.append(
                            f"line {ln} [{sid}] {col}: missing file '{f}'"
                        )

    if problems:
        raise WorkflowError(
            "[manifest] {} problem(s) in {}:\n  - ".format(len(problems), path)
            + "\n  - ".join(problems)
        )

    return df


MANIFEST = load_manifest(MANIFEST_PATH)

# sample_id -> row dict
SAMPLES = MANIFEST.set_index("sample_id", drop=False).to_dict(orient="index")
FAMILIES = sorted(MANIFEST["family"].unique())


def _ids(source=None, layout=None):
    sub = MANIFEST
    if source:
        sub = sub[sub["source"] == source]
    if layout:
        sub = sub[sub["layout"] == layout]
    return sub["sample_id"].tolist()


def _constraint(ids):
    """Regex matching exactly these ids, or nothing at all when the list is empty."""
    return "|".join(re.escape(i) for i in ids) if ids else r"(?!x)x"


SRA_PAIRED   = _ids("sra",   "paired")
SRA_SINGLE   = _ids("sra",   "single")
LOCAL_PAIRED = _ids("local", "paired")
LOCAL_SINGLE = _ids("local", "single")

# Global guard so {sample} / {family} cannot swallow path separators, and so
# the four staging rules below stay mutually exclusive.
wildcard_constraints:
    sample = _constraint(MANIFEST["sample_id"].tolist()),
    family = _constraint(FAMILIES),
    sra_id = r"[A-Za-z0-9]+",


# ---------------------------------------------------------------------
# Public helpers for use in downstream rules
# ---------------------------------------------------------------------
def sample_reads(sample_id):
    """Normalised read paths for one sample: 2 entries if paired, else 1."""
    row = SAMPLES[sample_id]
    base = f"{READS_DIR}/{row['family']}/{sample_id}/{sample_id}"
    if row["layout"] == "paired":
        return [f"{base}_1.fastq.gz", f"{base}_2.fastq.gz"]
    return [f"{base}.fastq.gz"]


def get_reads(wildcards):
    """
    Input function for any rule carrying a {sample} wildcard:

        rule fastp:
            input:  get_reads
            output: ...
    """
    return sample_reads(wildcards.sample)


def get_reads_named(wildcards):
    """
    Same as get_reads but as a dict, so paired rules can use {input.fq1} /
    {input.fq2}. Single-end samples get only 'fq1'.
    """
    files = sample_reads(wildcards.sample)
    d = {"fq1": files[0]}
    if len(files) == 2:
        d["fq2"] = files[1]
    return d


def family_samples(family):
    return MANIFEST.loc[MANIFEST["family"] == family, "sample_id"].tolist()


def get_family_read_targets(wildcards):
    """Every normalised read file belonging to one family."""
    targets = []
    for sid in family_samples(wildcards.family):
        targets += sample_reads(sid)
    return targets


# ---------------------------------------------------------------------
# SRA route
# ---------------------------------------------------------------------
rule fetch_and_dump_sra_paired:
    """Download one accession into the shared cache (paired-end)."""
    output:
        fq1 = _maybe_temp(f"{SRA_CACHE}/{{sra_id}}_1.fastq.gz"),
        fq2 = _maybe_temp(f"{SRA_CACHE}/{{sra_id}}_2.fastq.gz"),
    benchmark:
        "benchmarks/fetch_and_dump_sra_paired/{sra_id}.tsv"
    params:
        cache = SRA_CACHE,
    threads: 6
    # Capped: concurrent fasterq-dump jobs saturate the network before they
    # saturate CPU.
    resources: sra_fetch=1
    conda: get_env("ASSEMBLY_ENV")
    retries: 2
    shell:
        r"""
        # Failure counter persists across separate snakemake/sbatch
        # invocations, not just retries: 2's within-run attempts.
        FAIL_COUNT_FILE="$(pwd)/logs/sra_download_failures/{wildcards.sra_id}.count"
        MISMATCH_LOG="$(pwd)/logs/sra_layout_mismatches.tsv"
        mkdir -p "$(dirname "$FAIL_COUNT_FILE")"
        PRIOR_FAILS=0
        [ -f "$FAIL_COUNT_FILE" ] && PRIOR_FAILS=$(cat "$FAIL_COUNT_FILE")
        if [ "$PRIOR_FAILS" -ge 3 ]; then
            echo "ERROR: {wildcards.sra_id} has already failed $PRIOR_FAILS time(s) across previous attempts -- not retrying automatically. Investigate manually, then 'rm $FAIL_COUNT_FILE' to allow another attempt." >&2
            exit 1
        fi
        echo $((PRIOR_FAILS + 1)) > "$FAIL_COUNT_FILE"

        mkdir -p {params.cache}
        cd {params.cache}

        if aws s3 cp "s3://sra-pub-run-odp/sra/{wildcards.sra_id}/{wildcards.sra_id}" "{wildcards.sra_id}.sra" --no-sign-request 2>/dev/null; then
            SRA_SRC="{wildcards.sra_id}.sra"
        else
            prefetch --max-size 50G "{wildcards.sra_id}" || true
            SRA_SRC="{wildcards.sra_id}"
        fi
        # --split-3 (not --split-files) auto-detects layout: _1/_2 mates
        # for genuinely paired runs, one unsuffixed fastq for single-end.
        fasterq-dump --threads {threads} --split-3 --skip-technical \
                     "$SRA_SRC" -O .

        # This checkout's downstream is paired-only; a manifest row
        # claiming paired that's actually single-end fails loudly here
        # (logged once) rather than being faked with an empty mate file.
        if [ ! -f "{wildcards.sra_id}_2.fastq" ]; then
            [ -f "$MISMATCH_LOG" ] || echo -e "sra_id\tdetected_layout\tdetected_at" > "$MISMATCH_LOG"
            grep -q "^{wildcards.sra_id}	" "$MISMATCH_LOG" 2>/dev/null || \
                echo -e "{wildcards.sra_id}\tsingle\t$(date -Iseconds)" >> "$MISMATCH_LOG"
            echo "ERROR: {wildcards.sra_id} is manifest layout=paired but fasterq-dump --split-3 produced no _2.fastq -- this accession is actually single-end. Recorded in $MISMATCH_LOG. Update the manifest's layout column to 'single', then 'rm $FAIL_COUNT_FILE' to retry." >&2
            rm -rf "{wildcards.sra_id}"*.fastq "{wildcards.sra_id}" "{wildcards.sra_id}".sra "{wildcards.sra_id}".*.tmp
            exit 1
        fi

        # Concurrent, level 1,
        # half the threads each.
        HALF_THREADS=$(( {threads} / 2 )); [ "$HALF_THREADS" -lt 1 ] && HALF_THREADS=1
        pigz -1 -p "$HALF_THREADS" -c "{wildcards.sra_id}_1.fastq" > "{wildcards.sra_id}_1.fastq.gz" &
        PID1=$!
        pigz -1 -p "$HALF_THREADS" -c "{wildcards.sra_id}_2.fastq" > "{wildcards.sra_id}_2.fastq.gz" &
        PID2=$!
        wait "$PID1"
        wait "$PID2"
        rm -rf "{wildcards.sra_id}_1.fastq" "{wildcards.sra_id}_2.fastq" \
               "{wildcards.sra_id}" "{wildcards.sra_id}".sra "{wildcards.sra_id}".*.tmp

        rm -f "$FAIL_COUNT_FILE"
        """


rule fetch_and_dump_sra_single:
    """Download one accession into the shared cache (single-end / long read)."""
    output:
        fq = _maybe_temp(f"{SRA_CACHE}/{{sra_id}}.fastq.gz"),
    benchmark:
        "benchmarks/fetch_and_dump_sra_single/{sra_id}.tsv"
    params:
        cache = SRA_CACHE,
    threads: 6
    resources: sra_fetch=1
    conda: get_env("ASSEMBLY_ENV")
    retries: 2
    shell:
        r"""
        # See fetch_and_dump_sra_paired for the failure-counter rationale.
        FAIL_COUNT_FILE="$(pwd)/logs/sra_download_failures/{wildcards.sra_id}.count"
        mkdir -p "$(dirname "$FAIL_COUNT_FILE")"
        PRIOR_FAILS=0
        [ -f "$FAIL_COUNT_FILE" ] && PRIOR_FAILS=$(cat "$FAIL_COUNT_FILE")
        if [ "$PRIOR_FAILS" -ge 3 ]; then
            echo "ERROR: {wildcards.sra_id} has already failed $PRIOR_FAILS time(s) across previous attempts -- not retrying automatically. Investigate manually, then 'rm $FAIL_COUNT_FILE' to allow another attempt." >&2
            exit 1
        fi
        echo $((PRIOR_FAILS + 1)) > "$FAIL_COUNT_FILE"

        mkdir -p {params.cache}
        cd {params.cache}

        if aws s3 cp "s3://sra-pub-run-odp/sra/{wildcards.sra_id}/{wildcards.sra_id}" "{wildcards.sra_id}.sra" --no-sign-request 2>/dev/null; then
            SRA_SRC="{wildcards.sra_id}.sra"
        else
            prefetch --max-size 50G "{wildcards.sra_id}" || true
            SRA_SRC="{wildcards.sra_id}"
        fi
        fasterq-dump --threads {threads} --concatenate-reads --skip-technical \
                     "$SRA_SRC" -O .

        pigz -p {threads} -c "{wildcards.sra_id}.fastq" > "{wildcards.sra_id}.fastq.gz"
        rm -rf "{wildcards.sra_id}.fastq" "{wildcards.sra_id}" "{wildcards.sra_id}".sra \
               "{wildcards.sra_id}".*.tmp

        rm -f "$FAIL_COUNT_FILE"
        """


rule stage_sra_paired:
    input:
        fq1 = lambda wc: f"{SRA_CACHE}/{SAMPLES[wc.sample]['run_id']}_1.fastq.gz",
        fq2 = lambda wc: f"{SRA_CACHE}/{SAMPLES[wc.sample]['run_id']}_2.fastq.gz",
    output:
        fq1 = f"{READS_DIR}/{{family}}/{{sample}}/{{sample}}_1.fastq.gz",
        fq2 = f"{READS_DIR}/{{family}}/{{sample}}/{{sample}}_2.fastq.gz",
    benchmark:
        "benchmarks/stage_sra_paired/{family}_{sample}.tsv"
    wildcard_constraints:
        sample = _constraint(SRA_PAIRED),
    shell:
        """
        mkdir -p "$(dirname {output.fq1})"
        ln -sfn "$(readlink -f {input.fq1})" {output.fq1}
        ln -sfn "$(readlink -f {input.fq2})" {output.fq2}
        """


rule stage_sra_single:
    input:
        fq = lambda wc: f"{SRA_CACHE}/{SAMPLES[wc.sample]['run_id']}.fastq.gz",
    output:
        fq = f"{READS_DIR}/{{family}}/{{sample}}/{{sample}}.fastq.gz",
    benchmark:
        "benchmarks/stage_sra_single/{family}_{sample}.tsv"
    wildcard_constraints:
        sample = _constraint(SRA_SINGLE),
    shell:
        """
        mkdir -p "$(dirname {output.fq})"
        ln -sfn "$(readlink -f {input.fq})" {output.fq}
        """


# ---------------------------------------------------------------------
# Local route
# ---------------------------------------------------------------------
rule stage_local_paired:
    """Illumina paired-end fastq already on disk."""
    input:
        fq1 = lambda wc: _split_paths(SAMPLES[wc.sample]["fq1"]),
        fq2 = lambda wc: _split_paths(SAMPLES[wc.sample]["fq2"]),
    output:
        fq1 = f"{READS_DIR}/{{family}}/{{sample}}/{{sample}}_1.fastq.gz",
        fq2 = f"{READS_DIR}/{{family}}/{{sample}}/{{sample}}_2.fastq.gz",
    benchmark:
        "benchmarks/stage_local_paired/{family}_{sample}.tsv"
    wildcard_constraints:
        sample = _constraint(LOCAL_PAIRED),
    params:
        mode = LOCAL_LINK_MODE,
        script = STAGE_SCRIPT,
    shell:
        """
        bash {params.script} {params.mode} {output.fq1} {input.fq1}
        bash {params.script} {params.mode} {output.fq2} {input.fq2}
        """


rule stage_local_single:
    """ONT / PacBio / single-end fastq already on disk."""
    input:
        fq = lambda wc: _split_paths(SAMPLES[wc.sample]["fq1"]),
    output:
        fq = f"{READS_DIR}/{{family}}/{{sample}}/{{sample}}.fastq.gz",
    benchmark:
        "benchmarks/stage_local_single/{family}_{sample}.tsv"
    wildcard_constraints:
        sample = _constraint(LOCAL_SINGLE),
    params:
        mode = LOCAL_LINK_MODE,
        script = STAGE_SCRIPT,
    shell:
        "bash {params.script} {params.mode} {output.fq} {input.fq}"


# ---------------------------------------------------------------------
# Per-family collection point
# ---------------------------------------------------------------------
rule complete_family_reads:
    input:
        get_family_read_targets,
    output:
        done = f"{READS_DIR}/{{family}}/reads.done",
    benchmark:
        "benchmarks/complete_family_reads/{family}.tsv"
    run:
        rows = MANIFEST[MANIFEST["family"] == wildcards.family]
        with open(output.done, "w") as fh:
            fh.write(f"# family\t{wildcards.family}\n")
            fh.write(f"# n_samples\t{len(rows)}\n")
            fh.write("sample_id\tsource\tplatform\tlayout\trun_id\tbarcode\n")
            for _, r in rows.iterrows():
                fh.write("\t".join([
                    r["sample_id"], r["source"], r["platform"],
                    r["layout"], r["run_id"] or "-", r["barcode"] or "-",
                ]) + "\n")


# Back-compat: keeps `data/01.SRA_files/{family}/downloads.done` as a target.
rule complete_family_downloads:
    input:
        f"{READS_DIR}/{{family}}/reads.done",
    output:
        "data/01.SRA_files/{family}/downloads.done",
    benchmark:
        "benchmarks/complete_family_downloads/{family}.tsv"
    shell:
        """
        mkdir -p "$(dirname {output})"
        cp {input} {output}
        """
