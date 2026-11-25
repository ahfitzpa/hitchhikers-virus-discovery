def get_assembly_input(wildcards):
    return {
        "f1": f"data/02.fastp/{wildcards.family}/{wildcards.sample}/{wildcards.sample}_cleaned_1.fastq.gz",
        "f2": f"data/02.fastp/{wildcards.family}/{wildcards.sample}/{wildcards.sample}_cleaned_2.fastq.gz"
    }

# Highmem routing for assemblies that OOM on the shared partition, and the
# known-OOM exception list.
HIGHMEM_ATTEMPT = 3


def _highmem_samples():
    p = "configs/highmem_assembly_samples.txt"
    if not os.path.exists(p):
        return set()
    with open(p) as fh:
        return {ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")}


HIGHMEM_SAMPLES = _highmem_samples()


def _assembly_highmem(wildcards, attempt):
    return attempt >= HIGHMEM_ATTEMPT or wildcards.sample in HIGHMEM_SAMPLES


rule rnaviral_assembly:
    input:
        unpack(get_assembly_input)
    output:
        # temp(): contigs is consumed by filter_contigs; the BAM by
        # filter_contigs AND contig_coverage_stats, so Snakemake holds it until
        # the later of the two finishes. Coverage must be taken before the BAM
        # is reclaimed (for samples assembled before the coverage rule existed,
        # see remap_for_coverage).
        contigs = temp("data/03.rnaviralSPAdes/{family}/{sample}/contigs.fasta"),
        bam = temp("data/03.rnaviralSPAdes/{family}/{sample}/sorted.bam")
    benchmark:
        "benchmarks/rnaviral_assembly/{family}_{sample}.tsv"
    threads: 16
    # No mem_mb: this cluster's `shared` partition schedules on cores
    # only and rejects any --mem request outright; mem_gb below is
    # local-only bookkeeping for the free-memory probe in the shell block.
    resources:
        mem_gb=lambda wildcards, attempt: 256 if _assembly_highmem(wildcards, attempt) else 48,
        runtime=lambda wildcards, attempt: 960 if _assembly_highmem(wildcards, attempt) else 480,
        slurm_extra=lambda wildcards, attempt: ("--constraint=highmem"
                                                if _assembly_highmem(wildcards, attempt) else ""),
        # attempt_no is a resource, and the exception list is read in the
        # shell, because a params: block here would re-run every assembly
        attempt_no=lambda wildcards, attempt: attempt
    retries: 2
    # Not ASSEMBLY_ENV: bwa-mem3 and hostile=1.1.0's pinned
    # bowtie2 can't coexist in one conda env (zlib version conflict).
    conda: get_env("RNAVIRAL_ASSEMBLY_ENV")
    shell:
        """
        # A listed sample starts on highmem, so a retry would repeat it.
        HIGHMEM_LIST=configs/highmem_assembly_samples.txt
        if [ "{resources.attempt_no}" -gt 1 ] && [ -f "$HIGHMEM_LIST" ] \
           && grep -qx "{wildcards.sample}" "$HIGHMEM_LIST"; then
            echo "ERROR: {wildcards.sample} is in configs/highmem_assembly_samples.txt and its highmem attempt already failed." >&2
            echo "Not repeating an identical highmem run. Reduce its input instead (coverage normalisation or subsampling), then remove it from that list." >&2
            exit 1
        fi

        OUT_DIR=$(dirname "{output.contigs}")
        # Node-local disk, not BeeGFS.
        # SPAdes' own --tmp-dir only covers error-correction files, not the
        # per-k-mer scaffolding graphs that dominate its small-write volume, so
        # -o itself has to move.
        LOCAL_DIR="/tmp/rnaviral_assembly_{wildcards.family}_{wildcards.sample}_${{SLURM_JOB_ID:-$$}}"
        mkdir -p "$LOCAL_DIR"

        # Cleanup runs on EVERY exit, failure included. Before this, a failed
        # or timed-out SPAdes run left its whole -o dir (tens of GB) in the
        # node's /tmp, and later assemblies on that node died with "No space
        # left on device" (SRR16767818 on sonicmem5, 2026-11). The bwa-mem3
        # index (~1 GB/sample, 113 GB over the run) is only needed to build
        # sorted.bam and was never removed. sorted.bam.tmp.* are samtools sort
        # spill files a killed sort leaves behind.
        cleanup_assembly() {{
            rc=$?
            rm -rf "$LOCAL_DIR"
            rm -f "$OUT_DIR"/contigs.{{0123,amb,ann,bwt.2bit.64,pac}}
            if [ "$rc" -ne 0 ]; then
                rm -f "$OUT_DIR"/sorted.bam.tmp.*.bam
            fi
        }}
        trap cleanup_assembly EXIT
        # SLURM time limit / scancel send TERM: exit so the EXIT trap fires
        trap 'exit 143' TERM INT

        # Detect real free memory on whatever node this attempt lands on,
        # rather than trusting the static resources.mem_gb fallback blind.
        AVAIL_MEM_GB=$(( $(free -m | awk '/^Mem:/{{print $7}}') / 1024 ))
        if [ -z "$AVAIL_MEM_GB" ] || [ "$AVAIL_MEM_GB" -le 0 ]; then
            AVAIL_MEM_GB={resources.mem_gb}
        fi
        SPADES_MEM_GB=$(( AVAIL_MEM_GB / 2 ))
        if [ "$SPADES_MEM_GB" -lt 16 ]; then
            SPADES_MEM_GB=16
        fi
        if [ "$SPADES_MEM_GB" -gt {resources.mem_gb} ]; then
            SPADES_MEM_GB={resources.mem_gb}
        fi

        # Select only k-mers strictly shorter than the library's max read length
        # (prevents SPAdes exit status 67 on 75 bp libraries). Ported from
        # conifer_pest_virome (a42740f).
        #
        # rstrip() takes NO argument on purpose. Snakemake parses this shell
        # block as a Python string first, so a literal "\\r\\n" here would be
        # turned into a real CR+LF before the shell ever saw it, splitting the
        # Python string across lines ("SyntaxError: unterminated string
        # literal"). Bare rstrip() strips trailing whitespace including CR/LF
        # and needs no escaping.
        KMERS=$(python3 -c '
import gzip, sys
m = 0
with gzip.open(sys.argv[1], "rt") as fh:
    for i, line in enumerate(fh):
        if i >= 40000:
            break
        if i % 4 == 1:
            n = len(line.rstrip())
            if n > m:
                m = n
ks = [str(k) for k in (21, 33, 45, 55, 67, 77) if k < m]
print(",".join(ks) if ks else "auto")
' {input.f1})
        echo "[rnaviral_assembly] selected -k $KMERS" >&2

        rnaviralspades.py --pe1-1 {input.f1} --pe1-2 {input.f2} \
            -k "$KMERS" -o "$LOCAL_DIR" -t {threads} -m $SPADES_MEM_GB
        cp "$LOCAL_DIR/contigs.fasta" {output.contigs}
        # free node-local disk now rather than holding it through bwa-mem3
        rm -rf "$LOCAL_DIR"

        bwa-mem3 index -p "$OUT_DIR/contigs" {output.contigs}
        bwa-mem3 mem -t {threads} --bam=0 "$OUT_DIR/contigs" \
            {input.f1} {input.f2} | \
            samtools view -u -F 4 - | samtools sort -l 1 -m 1G -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule filter_contigs:
    input:
        contigs = "data/03.rnaviralSPAdes/{family}/{sample}/contigs.fasta",
        bam = "data/03.rnaviralSPAdes/{family}/{sample}/sorted.bam"
    output:
        filtered = "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta"
    benchmark:
        "benchmarks/filter_contigs/{family}_{sample}.tsv"
    # Threads are for seqkit grep only; samtools idxstats is single-threaded.
    # Not yet sacct-validated.
    threads: 4
    resources: mem_gb=4, runtime=30
    conda: get_env("ASSEMBLY_ENV")
    shell:
        """
        OUT_DIR=$(dirname "{input.contigs}")
        samtools idxstats {input.bam} | awk '$3==0 {{print $1}}' > "$OUT_DIR/zero_coverage_contigs.txt"
        seqkit grep -v -j {threads} -f "$OUT_DIR/zero_coverage_contigs.txt" {input.contigs} -o {output.filtered}
        # The .bai is NOT removed here any more: contig_coverage_stats also
        # reads this BAM (samtools idxstats needs the index) and removes the
        # .bai itself once it is done.
        """


# -------------------------------------------------------------------------
# Per-contig coverage (conifer_pest_virome capture spec item 4): mean, median
# and evenness. Ported from conifer, with one Serratus-only addition: the BAM
# can come from either of two places.
#
#   * sorted.bam -- the assembly's own BAM (reads mapped to every assembled
#     contig). Used for any sample assembled from now on: this rule is its
#     second declared consumer, so the temp() BAM survives until coverage is
#     taken.
#   * remap.bam  -- remap_for_coverage below. For samples assembled BEFORE
#     this rule existed, whose sorted.bam, contigs.fasta and reads were all
#     reclaimed long ago. Reads are fetched again and mapped to
#     filtered_contigs.fasta; nothing is reassembled.
#
# The choice is made at DAG-build time: an existing sorted.bam wins; otherwise
# a filtered_contigs.fasta that already exists means the assembly ran without
# a coverage consumer and its BAM is gone, so remap. A sample not yet
# assembled has neither and gets sorted.bam.
#
# Difference to expect between the two: remap.bam is mapped against the
# RETAINED contigs only (contigs.fasta no longer exists), so its stats cover
# retained contigs only, and reads that used to land on a dropped zero-coverage
# contig may now place elsewhere. Retained contigs' depth should barely move;
# the per-assembly "all contigs" totals in stage_counts/escalation_fractions
# cannot be recovered for these samples.
# -------------------------------------------------------------------------
def get_coverage_bam(wildcards):
    d = f"data/03.rnaviralSPAdes/{wildcards.family}/{wildcards.sample}"
    if os.path.exists(f"{d}/sorted.bam"):
        return f"{d}/sorted.bam"
    if os.path.exists(f"{d}/filtered_contigs.fasta"):
        return f"{d}/remap.bam"
    return f"{d}/sorted.bam"


# Reads for remap_for_coverage are fetched and cleaned HERE, not through
# fetch_and_dump_sra_paired -> stage_sra_paired -> fastp_cleaning. That chain
# also feeds rnaviral_assembly, and regenerating any file in it makes
# Snakemake treat the assembly's inputs as updated: a dry run of the first
# draft, which reused fastp_cleaning's output, planned all 96 novel assemblies
# and their whole taxonomy again (1353 jobs). This rule has no inputs, so
# nothing upstream of filtered_contigs.fasta is ever touched.
#
# The download and fastp calls are copied from fetch_and_dump_sra_paired
# (rules/01.input.smk) and fastp_cleaning (rules/02.fastp.smk), so the reads
# are cleaned exactly as they were for the assembly. Raw reads never leave
# this job; the cleaned pair is temp() and reclaimed after remap_for_coverage.
rule refetch_reads_for_remap:
    output:
        f1 = temp("data/02.fastp_remap/{family}/{sample}/{sample}_cleaned_1.fastq.gz"),
        f2 = temp("data/02.fastp_remap/{family}/{sample}/{sample}_cleaned_2.fastq.gz"),
        json = "data/02.fastp_remap/{family}/{sample}/{sample}_report.json"
    benchmark:
        "benchmarks/refetch_reads_for_remap/{family}_{sample}.tsv"
    params:
        sra_id = lambda wildcards: SAMPLES[wildcards.sample]["run_id"],
        min_len = 50,
        overlap = 30
    threads: 8
    resources: sra_fetch=1, runtime=480
    conda: get_env("ASSEMBLY_ENV")
    retries: 2
    shell:
        r"""
        OUT_DIR=$(dirname {output.f1})
        WORK="$OUT_DIR/.raw"
        trap 'rm -rf "$WORK"' EXIT
        trap 'exit 143' TERM INT
        rm -rf "$WORK"; mkdir -p "$WORK"
        cd "$WORK"
        if aws s3 cp "s3://sra-pub-run-odp/sra/{params.sra_id}/{params.sra_id}" "{params.sra_id}.sra" --no-sign-request 2>/dev/null; then
            SRA_SRC="{params.sra_id}.sra"
        else
            prefetch --max-size 50G "{params.sra_id}" || true
            SRA_SRC="{params.sra_id}"
        fi
        fasterq-dump --threads {threads} --split-3 --skip-technical "$SRA_SRC" -O .
        rm -rf "{params.sra_id}" "{params.sra_id}.sra"
        if [ ! -f "{params.sra_id}_2.fastq" ]; then
            echo "ERROR: {params.sra_id}: fasterq-dump --split-3 produced no _2.fastq (single-end?)" >&2
            exit 1
        fi
        cd - >/dev/null

        fastp --in1 "$WORK/{params.sra_id}_1.fastq" --in2 "$WORK/{params.sra_id}_2.fastq" \
            --out1 {output.f1} --out2 {output.f2} \
            --html /dev/null \
            --json {output.json} \
            --qualified_quality_phred 20 \
            --length_required {params.min_len} \
            --overlap_len_require {params.overlap} \
            --correction --overlap_diff_percent_limit 2 \
            -z 1 --thread {threads}
        """


rule remap_for_coverage:
    input:
        f1 = "data/02.fastp_remap/{family}/{sample}/{sample}_cleaned_1.fastq.gz",
        f2 = "data/02.fastp_remap/{family}/{sample}/{sample}_cleaned_2.fastq.gz",
        contigs = "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta"
    output:
        # temp(): contig_coverage_stats is the only consumer.
        bam = temp("data/03.rnaviralSPAdes/{family}/{sample}/remap.bam")
    benchmark:
        "benchmarks/remap_for_coverage/{family}_{sample}.tsv"
    threads: 16
    resources: mem_gb=32, runtime=480
    # Same env and the same bwa-mem3 / samtools calls as rnaviral_assembly, so
    # the two BAM sources differ only in the reference.
    conda: get_env("RNAVIRAL_ASSEMBLY_ENV")
    shell:
        """
        OUT_DIR=$(dirname {output.bam})
        IDX="$OUT_DIR/remap_idx"
        cleanup_remap() {{
            rc=$?
            rm -f "$IDX".{{0123,amb,ann,bwt.2bit.64,pac}}
            if [ "$rc" -ne 0 ]; then
                rm -f {output.bam}.tmp.*.bam {output.bam}.bai
            fi
        }}
        trap cleanup_remap EXIT
        trap 'exit 143' TERM INT

        bwa-mem3 index -p "$IDX" {input.contigs}
        bwa-mem3 mem -t {threads} --bam=0 "$IDX" \
            {input.f1} {input.f2} | \
            samtools view -u -F 4 - | samtools sort -l 1 -m 1G -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule contig_coverage_stats:
    input:
        bam = get_coverage_bam,
        filtered = "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta"
    output:
        tsv = "data/03.rnaviralSPAdes/{family}/{sample}/coverage_stats.tsv",
        idxstats = "data/03.rnaviralSPAdes/{family}/{sample}/idxstats.tsv"
    benchmark:
        "benchmarks/contig_coverage_stats/{family}_{sample}.tsv"
    threads: 2
    resources: mem_gb=8, runtime=60
    conda: get_env("ASSEMBLY_ENV")
    shell:
        # -a: report zero-depth positions too, so breadth and evenness are
        # computed over the full contig length rather than covered bases only.
        # -Q 0 -q 0: no MAPQ/base-quality gate, matching the unfiltered
        # idxstats counts filter_contigs uses.
        #
        # No --max-depth: samtools >=1.13 dropped the depth cap flag and errors
        # on it; unlimited depth is the default.
        """
        samtools idxstats {input.bam} > {output.idxstats}
        samtools depth -a -Q 0 -q 0 {input.bam} \
            | python scripts/contig_coverage_stats.py \
                --idxstats {output.idxstats} \
                --depth - \
                --retained {input.filtered} \
                --family {wildcards.family} \
                --sample {wildcards.sample} \
                --out {output.tsv}
        # The BAM is temp() and reclaimed after this rule; its index is not a
        # declared output, so remove it here.
        rm -f {input.bam}.bai
        """


# Annotation only since the three-tool taxonomy switch: its score is joined
# onto every sample's taxonomy_merged.tsv but no longer gates what the
# taxonomy tools see. The COBRA stage that used to follow it was removed.
rule run_virsorter2:
    # runtime: was 360 (6h), which 142 samples clipped in the 2026-11-02 run and
    # which then aborted the workflow at 86%. Benchmarks over 84 completed runs:
    # median 1.97h, p90 4.86h, max 5.71h -- nothing exceeded 6h, so the limit was
    # marginal rather than the tool hanging. 900 leaves real headroom.
    conda: get_env("TAX_QC_ENV")
    threads: 6
    resources: runtime=900
    retries: 2
    input:
        fasta = "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta",
        db = os.path.join(config["databases"]["virsorter2_db"], ".done")
    output:
        tsv = "data/03.rnaviralSPAdes/{family}/{sample}/virsorter2/final-viral-score.tsv",
        fasta = "data/03.rnaviralSPAdes/{family}/{sample}/virsorter2/final-viral-combined.fa",
        dir = directory("data/03.rnaviralSPAdes/{family}/{sample}/virsorter2")
    benchmark:
        "benchmarks/run_virsorter2/{family}_{sample}.tsv"
    shell:
        """
        # Best-effort: both outputs are always created, empty on failure, and
        # failures are logged.
        mkdir -p {output.dir} logs
        # Run under `timeout`, not bare. SLURM's wall-clock kill terminates this
        # shell outright, so the handler below never ran and a TIMEOUT failed the
        # rule hard -- which is what held this workflow in stage 1 for 1 d 21 h
        # across 17 timeouts on six degenerate assemblies (1.7-4.7M contigs
        # against a 147k median). Owning the kill ourselves means the command
        # returns 124, the handler runs, empty outputs are written and the sample
        # degrades to "no score" exactly as a crash already did. The 600 s margin
        # leaves room for the output guarantees and the working-file purge below.
        VS2_LIMIT=$(( {resources.runtime} * 60 ))
        # Margin is 600 s, except for short allocations where a flat 600 s would
        # eat most or all of the budget (and go negative under ~10 min), leaving
        # no room for the handler -- there, take a tenth.
        if [ "$VS2_LIMIT" -gt 6000 ]; then VS2_MARGIN=600; else VS2_MARGIN=$(( VS2_LIMIT / 10 )); fi
        VS2_BUDGET=$(( VS2_LIMIT - VS2_MARGIN ))
        VS2_RC=0
        timeout -k 60 "${{VS2_BUDGET}}s" \
            virsorter run -i {input.fasta} -w {output.dir} \
                --db-dir {config[databases][virsorter2_db]} -j {threads} || VS2_RC=$?
        if [ "$VS2_RC" -ne 0 ]; then
            # 124 = timeout sent SIGTERM; 137 = it escalated to SIGKILL.
            if [ "$VS2_RC" -eq 124 ] || [ "$VS2_RC" -eq 137 ]; then
                VS2_REASON="TIMEOUT after ${{VS2_BUDGET}}s"
            else
                VS2_REASON="FAILED rc=$VS2_RC"
            fi
            echo "[virsorter2] $VS2_REASON for {wildcards.family}/{wildcards.sample}; continuing without its score" >&2
            printf '%s\t%s\t%s\t%s\n' "{wildcards.family}" "{wildcards.sample}" "$(date -Is)" "$VS2_REASON" \
                >> logs/virsorter2_failures.tsv
        fi
        # guarantee the declared outputs exist, including after a partial run
        [ -s {output.tsv} ] || : > {output.tsv}
        [ -f {output.fasta} ] || : > {output.fasta}
        # Nothing reads VirSorter2's working files (iter-0/ alone was ~1 GB a
        # sample, 210 GB over the run). Keep final-viral-*, config.yaml, log/.
        find {output.dir} -mindepth 1 -maxdepth 1 \
            ! -name 'final-viral-*' ! -name 'config.yaml' ! -name 'log' \
            ! -name '.snakemake_timestamp' -exec rm -rf {{}} +
        """

rule seqkit_stats:
    # Was on raw reads (data/01.reads); moved to the rnaviralSPAdes assembly
    # output per-sample instead. filtered_contigs.fasta, not
    # contigs.fasta: the latter is temp-wrapped and deleted once filter_contigs
    # finishes.
    input:
        "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta"
    conda: get_env("ASSEMBLY_ENV")
    output:
        "data/03.assembly_stats/{family}/{sample}/{sample}_stats.txt"
    benchmark:
        "benchmarks/seqkit_stats/{family}_{sample}.tsv"
    shell:
        """
        mkdir -p $(dirname {output})
        seqkit stats {input} > {output}
        """
