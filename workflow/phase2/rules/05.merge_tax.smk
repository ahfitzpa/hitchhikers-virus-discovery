# Per-sample union of the three taxonomy tools (rules/04.taxonomy.smk),
# then CheckV + VirSorter2 annotation, then the cross-sample compile.
# Column layout is documented in scripts/merge_taxonomy_calls.py.

rule collect_viral_positives:
    # Union, not intersection: a contig called viral by ANY of Kraken2,
    # geNomad or RdRpCATCH is kept. detected_by/n_tools record which.
    conda: get_env("TAX_CORE_ENV")
    # Seconds of work, but deliberately not a localrule. A
    # short runtime backfills queueing behind the 180 min default.
    threads: 1
    resources: mem_gb=4, runtime=30
    input:
        contigs=get_taxonomy_contigs,
        kraken2="data/04.taxonomy/{family}/{sample}/kraken2/kraken2_output.tsv",
        genomad="data/04.taxonomy/{family}/{sample}/genomad/filtered_contigs_summary/filtered_contigs_virus_summary.tsv",
        rdrpcatch="data/04.taxonomy/{family}/{sample}/rdrpcatch/rdrp_annotated.tsv",
        ktaxonomy=os.path.join(config["databases"]["kraken2_db"], "ktaxonomy.tsv")
    output:
        calls="data/05.final_selection/{family}/{sample}/tool_calls.tsv",
        fasta="data/05.final_selection/{family}/{sample}/viral_positives.fasta"
    benchmark:
        "benchmarks/collect_viral_positives/{family}_{sample}.tsv"
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/merge_taxonomy_calls.py collect \\
            --contigs {input.contigs} --kraken2 {input.kraken2} --ktaxonomy {input.ktaxonomy} \\
            --genomad {input.genomad} --rdrpcatch {input.rdrpcatch} \\
            --out-tsv {output.calls} --out-fasta {output.fasta}
        """

rule annotate_viral_positives:
    # CheckV quality and VirSorter2 score are joined on as columns only --
    # neither removes a contig the three taxonomy tools called positive.
    conda: get_env("TAX_CORE_ENV")
    threads: 1
    resources: mem_gb=4, runtime=30
    input:
        calls="data/05.final_selection/{family}/{sample}/tool_calls.tsv",
        checkv="data/04.taxonomy/{family}/{sample}/checkv/quality_summary.tsv"
    params:
        # params, not input: as an input a VirSorter2 failure blocked this
        # sample's merge. Ordering comes from stage 1 of submit_phase2.sh.
        virsorter2="data/03.rnaviralSPAdes/{family}/{sample}/virsorter2/final-viral-score.tsv"
    output:
        merged="data/05.final_selection/{family}/{sample}/taxonomy_merged.tsv"
    benchmark:
        "benchmarks/annotate_viral_positives/{family}_{sample}.tsv"
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/merge_taxonomy_calls.py annotate \\
            --calls {input.calls} --checkv {input.checkv} --virsorter2 {params.virsorter2} \\
            --out {output.merged}
        """

def get_all_final_files(wildcards):
    return [
        f"data/05.final_selection/{row.family}/{row.sample_id}/taxonomy_merged.tsv"
        for _, row in MANIFEST.iterrows()
    ]


# Partial aggregation: with PHASE2_ALLOW_PARTIAL_MASTER set (stage 3 of
# submit_phase2.sh), aggregate the merges that exist and record the rest in
# excluded_samples.tsv, instead of one unassemblable sample blocking the master
# database and MetaCerberus behind it. Default off: in stage 2 a missing merge
# is still a failure.
def _allow_partial_master():
    """Env var, not --config: this Snakefile discards --config entirely."""
    if str(os.environ.get("PHASE2_ALLOW_PARTIAL_MASTER", "")).lower() in ("1", "true", "yes"):
        return True
    return bool(config.get("allow_partial_master", False))


def get_master_inputs(wildcards):
    want = get_all_final_files(wildcards)
    if not _allow_partial_master():
        return want
    return [f for f in want if os.path.exists(f)
            and os.path.exists(f.replace("taxonomy_merged.tsv", "viral_positives.fasta"))]

# Per-tool columns (kraken2_*/genomad_*/rdrpcatch_*) plus consensus_family;
# keyed on (sample, contig_id), since SPAdes NODE_* names repeat across samples.
rule build_master_database:
    input:
        tables = get_master_inputs,
        fastas = lambda wildcards: [f.replace("taxonomy_merged.tsv", "viral_positives.fasta") for f in get_master_inputs(wildcards)],
        meta_dict = config["metadata"]["shortlist"],
        # the script runs outside the Snakefile, so it reads the manifest
        # itself rather than relying on the MANIFEST global
        manifest = MANIFEST_PATH
    output:
        master_fa = "results/master_viral_database.fasta",
        master_csv = "results/taxonomy/master_metadata.csv",
        # One row per manifest sample left out; always written.
        excluded = "results/taxonomy/excluded_samples.tsv"
    benchmark:
        "benchmarks/build_master_database/build.tsv"
    # TAX_CORE_ENV, not ASSEMBLY_ENV: the script imports pandas and Bio, and
    # env_phase2_assembly.yaml ships neither (awscli/sra-tools/fastp/seqkit/
    # rsync/samtools only). Moving this rule from `run:` to `script:` fixed
    # only half the bug -- in the assembly env it still dies on
    # ModuleNotFoundError. tax_core.yaml already provides pandas and biopython.
    # Fixed 2026-11-17, before this rule had ever run here.
    conda: get_env("TAX_CORE_ENV")
    # A `script:`, not a `run:`. Snakemake ignores `conda:` for `run:` blocks and
    # executes them in its own interpreter, so `from Bio import SeqIO` raised
    # ModuleNotFoundError and this rule could never succeed.
    script:
        "../scripts/build_master_database.py"


rule taxonomy_cleanup:
    # Purges tool scratch only -- never a file another rule declares as
    # output (see the output-tracking note atop rules/04.taxonomy.smk).
    input:
        final="data/05.final_selection/{family}/{sample}/taxonomy_merged.tsv"
    output:
        "data/04.taxonomy/{family}/{sample}/.purged"
    benchmark:
        "benchmarks/taxonomy_cleanup/{family}_{sample}.tsv"
    shell:
        """
        BASE="data/04.taxonomy/{wildcards.family}/{wildcards.sample}"
        rm -rf "$BASE/rdrpcatch/scan"
        rm -rf "$BASE/checkv/tmp"
        touch {output}
        """
