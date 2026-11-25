# Item 7: ICTV phylogenetic comparison per virus family -- SCAFFOLD ONLY.
#
# The method was explicitly deferred to a later prompt, so this file defines the
# interface and nothing more:
#
#   ictv_prepare          groups retained viral contigs by assigned_family and
#                         writes one self-contained bundle per family
#                         (sequences, RdRp domains, metadata, a TODO listing
#                         what the method still needs). Fully implemented.
#
#   ictv_family_comparison  the placeholder for the comparison itself. It
#                         declares the inputs and outputs the method will use
#                         and FAILS LOUDLY if requested, rather than emitting an
#                         empty result that could be mistaken for a real one.
#                         It is deliberately NOT in rule all.
#
# Family bundles are discovered at runtime (families depend on what the taxonomy
# tools call), so ictv_prepare writes a directory() plus a families.tsv index
# rather than per-family declared outputs. When the method is specified, convert
# ictv_prepare to a checkpoint and fan out over families.tsv.

ICTV_DIR = "results/ictv"


rule ictv_prepare:
    input:
        candidates = REPORT_DIR + "/candidate_contigs.tsv",
        master = "results/master_viral_database.fasta",
        rdrp = get_all_rdrp_domains
    output:
        families = ICTV_DIR + "/families.tsv",
        bundles = directory(ICTV_DIR + "/families")
    benchmark:
        "benchmarks/ictv_prepare/prepare.tsv"
    threads: 1
    resources: mem_gb=8, runtime=120
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/ictv_prepare.py \
            --candidates {input.candidates} \
            --master-fasta {input.master} \
            --rdrp-faa {input.rdrp} \
            --out-dir {output.bundles} \
            --families-tsv {output.families}
        """


# Placeholder. Not reachable from rule all -- requesting it is an error until
# the method is specified.
rule ictv_family_comparison:
    input:
        bundle = ICTV_DIR + "/families/{virus_family}/.prepared",
        contigs = ICTV_DIR + "/families/{virus_family}/contigs.fasta",
        rdrp = ICTV_DIR + "/families/{virus_family}/rdrp_domains.faa",
        # To be added when specified: the ICTV reference set for this family.
    output:
        tree = ICTV_DIR + "/families/{virus_family}/tree.newick",
        distances = ICTV_DIR + "/families/{virus_family}/ictv_distances.tsv",
        summary = ICTV_DIR + "/families/{virus_family}/ictv_summary.tsv"
    benchmark:
        "benchmarks/ictv_family_comparison/{virus_family}.tsv"
    threads: 8
    resources: mem_gb=16, runtime=480
    shell:
        """
        echo "ictv_family_comparison is a scaffold: no method specified yet." >&2
        echo "Needed: ICTV release + reference set, marker (genome / RdRp /" >&2
        echo "named protein), alignment + trimming, tree method + model," >&2
        echo "demarcation criterion, and the reported distance WITH UNITS." >&2
        echo "See {input.bundle%/.prepared}/REFERENCES.todo" >&2
        exit 1
        """
