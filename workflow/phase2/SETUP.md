# Phase 2 database setup

Phase 2 needs six reference databases before it can run `rule all`. This
is a one-time provisioning step per checkout, or per shared location if
several people/checkouts share one install (see the note at the bottom).
On UCD HPC Sonic they are already installed under `/path/to/databases/`
and `configs/paths.yaml` points there, so there is nothing to do.

| Database | Used by | Approx. size | Approx. time |
|---|---|---|---|
| Kraken2 viral (`k2_viral`) | `run_kraken2` | ~2GB | minutes |
| geNomad (v1.9) | `run_genomad` | ~1.5GB | minutes |
| RdRpCATCH | `run_rdrpcatch` | ~4GB | ~5-10 min |
| CheckV | `run_checkv` | ~6GB | ~5-10 min |
| VirSorter2 | `run_virsorter2` | ~12GB | ~10-20 min |
| MetaCerberus HMMs (KOFam_all, VOG, PHROG) | `run_metacerberus` | ~10-20GB | ~30 min |

The MMseqs2 `nt` database (~14TB on disk) is no longer used by phase2's
taxonomy (replaced by Kraken2 + geNomad + RdRpCATCH on 2026-11-02), so it is
not provisioned here.

## Prerequisites

- Snakemake 8.x (`snakemake --version`) and conda/mamba on your `PATH`.
- **Outbound internet access on whatever host actually runs this.** Every
  one of these downloads talks to a public server (NCBI, the VirSorter2/
  CheckV/RdRpCATCH DB mirrors, GitHub). On HPC clusters, compute nodes
  often can't reach the internet even though login nodes can — confirmed
  true for login nodes on this project's own cluster, unconfirmed for its
  compute nodes. If in doubt, run this from a login/head node, or
  whichever host you already know has a normal internet connection.
- **~50GB free disk.** Point it at scratch/bulk storage, not a small
  home-directory quota.

## Quick start

```bash
cd workflow/phase2
./setup_databases.sh /path/to/database/root
```

This points `configs/paths.yaml`'s `databases:` block at
`/path/to/database/root/<name>` for each of the six databases, then runs
each `setup_*` rule in turn (cheapest first, so a missing tool or bad
path fails in minutes).

The whole run can take an hour or more, so don't run it in a foreground
terminal you might lose. Either:

```bash
nohup ./setup_databases.sh /path/to/database/root > setup.log 2>&1 &
```

or, on a SLURM cluster where the login node has internet access but you'd
rather not tie up your own session:

```bash
sbatch --wrap="./setup_databases.sh /path/to/database/root" \
    --job-name=db_setup --cpus-per-task=8 --time=08:00:00 \
    --output=setup_%j.log --error=setup_%j.err
```

(Only do this if you've confirmed the partition's compute nodes actually
have internet access — otherwise the downloads will hang or fail there
even though the job itself runs fine.)

## Resuming after an interruption

Every `setup_*` rule is idempotent: it checks for its own `.done` flag (or
the specific files it needs, e.g. `hash.k2d` for Kraken2)
before doing any work, so re-running `setup_databases.sh` after a
network blip or a killed job just picks up wherever it left off — it
won't re-download anything already present.

## Manual / one-database-at-a-time alternative

If you'd rather not run the wrapper script (e.g. you already have some of
these databases elsewhere and just want the rest), edit
`configs/paths.yaml`'s `databases:` block yourself, then run whichever
rules you need directly:

```bash
snakemake -s Snakefile_phase2 --until setup_genomad \
    --cores 8 --use-conda --conda-prefix .snakemake/conda --conda-frontend conda
```

(swap `setup_genomad` for `setup_kraken2`, `setup_rdrpcatch`,
`setup_checkv`, `setup_virsorter2` or `setup_metacerberus` as needed).

## Sharing one install across multiple checkouts

If you're setting this up for a lab/group where several people will each
clone this repo (or the `deer_virome`-style forks that reuse it), install
the databases **once** at a shared path and point every checkout's
`configs/paths.yaml databases:` block at that same location — see the
main `README.md`'s "Shared databases" section for why this is the
intended pattern here, not a one-off optimization.
