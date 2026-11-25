import os

# Tier the viral positives into novel-candidate confidence classes.
#
# Thin wrapper around scripts/tier_novel_candidates.sh, which does the work in
# bash + awk and can also be run by hand against any --root while the pipeline
# is still going. The rule exists so the table is produced and dated alongside
# every other phase2 result rather than reconstructed ad hoc afterwards.
#
# This is a REPORTING step, not a gate: nothing downstream consumes its output,
# and it never removes a contig. Tiers and flags are advisory, and every row
# carries the reasons behind its tier so the call can be audited.
#
# It is NOT dereplication. Identical viruses recovered from several libraries
# are several rows here; --exact-dup-check catches only byte-identical contigs
# (re-submitted runs). Cluster the markers with mmseqs2/cd-hit before counting
# distinct viruses, or every total is inflated -- on the Hepeviridae set twelve
# ERR libraries carried RdRp blocks identical at the amino acid level while
# their contigs differed in length, so none of them were caught here.


def _novelty_genome_window(family):
    """(min, max) plausible complete-genome length for a family."""
    nov = config.get("novelty", {})
    window = (nov.get("genome_len", {}) or {}).get(family)
    if window and len(window) == 2:
        return int(window[0]), int(window[1])
    return int(nov.get("default_genome_min", 3000)), int(nov.get("default_genome_max", 20000))


def get_novelty_inputs(wildcards):
    """Every merge for this family, or just those that exist in partial mode.

    Mirrors build_master_database's tolerance (rules/05.merge_tax.smk) so that
    one unassemblable sample cannot withhold the table for the other 190.
    """
    want = [
        f"data/05.final_selection/{wildcards.family}/{s}/taxonomy_merged.tsv"
        for s in family_samples(wildcards.family)
    ]
    if not _allow_partial_master():
        return want
    return [f for f in want if os.path.exists(f)]


rule novel_candidate_tiers:
    input:
        tables = get_novelty_inputs,
        script = "scripts/tier_novel_candidates.sh"
    output:
        tsv = "results/reporting/novel_candidates_{family}.tsv"
    params:
        min_len = lambda wc: int(config.get("novelty", {}).get("min_contig_len", 1000)),
        gmin = lambda wc: _novelty_genome_window(wc.family)[0],
        gmax = lambda wc: _novelty_genome_window(wc.family)[1],
        dup = lambda wc: ("--exact-dup-check"
                          if config.get("novelty", {}).get("exact_dup_check", True) else ""),
    benchmark:
        "benchmarks/novel_candidate_tiers/{family}.tsv"
    threads: 1
    # Pure bash + awk over the per-sample TSVs: ~12 s for 191 Hepeviridae
    # samples, and no conda env (hence no `conda:` directive -- adding one
    # would build an environment this rule has no use for).
    resources: mem_gb=4, runtime=60
    log:
        "logs/novel_candidate_tiers/{family}.log"
    shell:
        # The script writes its tier/flag summary to stderr and the table to
        # --out, so the log keeps the summary next to the result.
        """
        bash {input.script} \
            --root . \
            --family {wildcards.family} \
            --min-len {params.min_len} \
            --genome-min {params.gmin} \
            --genome-max {params.gmax} \
            {params.dup} \
            --out {output.tsv} 2> {log}
        cat {log} >&2
        """
