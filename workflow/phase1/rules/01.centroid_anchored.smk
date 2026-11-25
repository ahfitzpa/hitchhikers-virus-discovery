# Slices must stay contiguous, not
# interleaved -- combine_fasta concatenates chunk_001.chunk_NNN in order and
# dedup keeps the first occurrence, so the original run_id order has to survive
# chunking.
N_DOWNLOAD_CHUNKS = 24
DOWNLOAD_CHUNK_IDS = [f"{i:03d}" for i in range(1, N_DOWNLOAD_CHUNKS + 1)]

def _get_run_ids(mode):
    checkpoint_output = checkpoints.centroid_step1_query.get(mode=mode).output["runs"]
    import pandas as pd
    return pd.read_csv(checkpoint_output)['run_id'].tolist()

def _chunk_slice(run_ids, n_chunks, chunk_num):
    n = len(run_ids)
    base, rem = divmod(n, n_chunks)
    start = chunk_num * base + min(chunk_num, rem)
    end = start + base + (1 if chunk_num < rem else 0)
    return run_ids[start:end]

def get_download_chunk_run_ids(wildcards):
    run_ids = _get_run_ids(wildcards.mode)
    return _chunk_slice(run_ids, N_DOWNLOAD_CHUNKS, int(wildcards.chunk) - 1)

rule centroid_download_chunk:
    input:
        runs=lambda wc: checkpoints.centroid_step1_query.get(mode=wc.mode).output["runs"]
    output: temp("data/01.sOTUs_centroid/{mode}/download_chunks/chunk_{chunk}.fa")
    log: "logs/{mode}/download_chunk_{chunk}.log"
    conda: "../" + config["tools"]["DISCOVERY_ENV"]
    params:
        run_ids=get_download_chunk_run_ids
    shell:
        """
        set -o pipefail
        : > {log}
        : > {output}
        REAL_ERROR=0
        ERR_LOG="{resources.tmpdir}/centroid_download_{wildcards.mode}_{wildcards.chunk}.err"
        for RUN_ID in {params.run_ids}; do
            S3="s3://logan-pub/u/${{RUN_ID}}/${{RUN_ID}}.unitigs.fa.zst"
            : > "$ERR_LOG"
            if aws s3 cp "$S3" - --no-sign-request 2>"$ERR_LOG" | zstd -d >> {output} 2>>{log}; then
                :
            elif grep -qE 'NoSuchKey|404' "$ERR_LOG"; then
                # Genuinely not in Logan -- expected, not an error. Write
                # nothing for this run_id and move on, same as before.
                echo "NOTE: $RUN_ID not found in Logan (NoSuchKey/404)" >> {log}
            else
                # Any other failure (throttling, network blip, a real S3 error)
                # must not be silently treated the same as "not in Logan" /
                # plan step 14: the old ls-then-touch pattern this replaced
                # made exactly this mistake and turned transient errors into
                # silently empty, accuracy-wrong output.
                echo "ERROR: $RUN_ID download failed for a reason other than NoSuchKey/404:" >> {log}
                cat "$ERR_LOG" >> {log}
                REAL_ERROR=1
            fi
        done
        rm -f "$ERR_LOG"
        if [ "$REAL_ERROR" -eq 1 ]; then
            echo "ERROR: at least one run_id in this chunk failed for a reason other than NoSuchKey/404 -- see {log}" >&2
            exit 1
        fi
        """

rule combine_fasta:
    input:
        expand(
            "data/01.sOTUs_centroid/{{mode}}/download_chunks/chunk_{chunk}.fa",
            chunk=DOWNLOAD_CHUNK_IDS
        )
    output: "data/01.sOTUs_centroid/{mode}/combined.fasta"
    shell:
        "cat {input} > {output} || touch {output}"

# Build the MMseqs query DB (deduplicated, clustered Logan unitigs) in three
# stages so chunk clustering can run as independent SLURM jobs. for why this is
# split this way and the fix history each rule below carries forward.
N_CENTROID_CHUNKS = 16
CENTROID_CHUNK_IDS = [f"{i:03d}" for i in range(1, N_CENTROID_CHUNKS + 1)]

rule centroid_dedup_and_split:
    input: "data/01.sOTUs_centroid/{mode}/combined.fasta"
    output:
        chunks = expand(
            "data/01.sOTUs_centroid/{mode}/dedup/dedup.part_{chunk}.fasta",
            chunk=CENTROID_CHUNK_IDS, allow_missing=True
        ),
        exact_dups = "data/01.sOTUs_centroid/{mode}/dedup/exact_dups.tsv"
    threads: 16
    resources:
        local_mem_mb=128000,
        runtime=2880
    conda: "../" + config["tools"]["DISCOVERY_ENV"]
    shell:
        """
        set -euo pipefail
        export LC_ALL=C

        OUT_DIR="data/01.sOTUs_centroid/{wildcards.mode}/dedup"
        BUILD_DIR="{resources.tmpdir}/centroid_dedup_{wildcards.mode}"
        rm -rf "$BUILD_DIR" "$OUT_DIR"
        mkdir -p "$BUILD_DIR" "$OUT_DIR"

        # Size the sort buffer off real free memory, not the static resource
        # declaration.
        AVAIL_MEM_MB=$(free -m | awk '/^Mem:/{{print $7}}')
        if [ -z "$AVAIL_MEM_MB" ] || [ "$AVAIL_MEM_MB" -le 0 ]; then
            AVAIL_MEM_MB={resources.local_mem_mb}
        fi
        SORT_BUF_MB=$(( AVAIL_MEM_MB / 2 ))
        SORT_BUF_CAP_MB=64000
        if [ "$SORT_BUF_MB" -gt "$SORT_BUF_CAP_MB" ]; then
            SORT_BUF_MB=$SORT_BUF_CAP_MB
        fi

        # idx (original order) / id / sequence / hash, one row per record.
        seqkit fx2tab --only-id --seq-hash -Q -j {threads} {input} \
            | awk -F'\t' '{{print NR"\t"$0}}' > "$BUILD_DIR/records.tsv"

        sort -t $'\t' -k4,4 -k1,1n -S "${{SORT_BUF_MB}}M" -T "$BUILD_DIR" \
            --parallel={threads} --compress-program=pigz \
            "$BUILD_DIR/records.tsv" > "$BUILD_DIR/records.sorted.tsv"

        touch "$BUILD_DIR/dedup.fasta" "{output.exact_dups}"
        awk -F'\t' -v fasta_out="$BUILD_DIR/dedup.fasta" -v dups_out="{output.exact_dups}" '
            function flush() {{
                if (n == 0) return
                print ">"ids[1] > fasta_out
                print seq > fasta_out
                for (i = 2; i <= n; i++) print ids[1]"\t"ids[i] > dups_out
            }}
            $4 != hash {{ flush(); hash = $4; seq = $3; n = 0 }}
            {{ n++; ids[n] = $2 }}
            END {{ flush() }}
        ' "$BUILD_DIR/records.sorted.tsv"

        rm -f "$BUILD_DIR/records.tsv" "$BUILD_DIR/records.sorted.tsv"

        seqkit split2 -p {N_CENTROID_CHUNKS} -O "$OUT_DIR" -j {threads} "$BUILD_DIR/dedup.fasta"
        rm -rf "$BUILD_DIR"
        """

rule centroid_cluster_chunk:
    input: "data/01.sOTUs_centroid/{mode}/dedup/dedup.part_{chunk}.fasta"
    output:
        repdb = "data/01.sOTUs_centroid/{mode}/chunks/{chunk}.repDB.dbtype",
        repdb_h = "data/01.sOTUs_centroid/{mode}/chunks/{chunk}.repDB_h.dbtype",
        clusters = "data/01.sOTUs_centroid/{mode}/chunks/{chunk}.clusters.tsv"
    wildcard_constraints:
        chunk = "[0-9]{3}"
    threads: 8
    resources:
        local_mem_mb=64000,
        runtime=1440
    conda: "../" + config["tools"]["DISCOVERY_ENV"]
    params:
        clust = lambda wildcards: config["querydb_clustering"]
    shell:
        """
        set -euo pipefail
        export LC_ALL=C

        verify_mmseqs_db() {{
            local db="$1" tries=0
            while [ ! -s "$db.dbtype" ] || [ ! -s "$db.index" ]; do
                tries=$((tries + 1))
                if [ "$tries" -ge 10 ]; then
                    echo "ERROR: $db is missing or incomplete after ${{tries}} checks over ~20s." >&2
                    return 1
                fi
                sleep 2
            done
        }}

        OUT_DIR="data/01.sOTUs_centroid/{wildcards.mode}/chunks"
        mkdir -p "$OUT_DIR"
        REPDB="$OUT_DIR/{wildcards.chunk}.repDB"
        BUILD_DIR="{resources.tmpdir}/centroid_chunk_{wildcards.mode}_{wildcards.chunk}"
        rm -rf "$BUILD_DIR"
        mkdir -p "$BUILD_DIR"
        CHUNK_DB="$BUILD_DIR/chunk.db"
        CHUNK_CLUSTERDB="$BUILD_DIR/chunk.clusterDB"
        CHUNK_TMP="$BUILD_DIR/chunk.tmp"

        mmseqs createdb {input} "$CHUNK_DB" \
            --createdb-mode 1 --shuffle 0 --write-lookup 0 \
            --threads {threads} -v 1

        mmseqs linclust "$CHUNK_DB" "$CHUNK_CLUSTERDB" "$CHUNK_TMP" \
            --threads {threads} \
            --min-seq-id {params.clust[min_seq_id]} \
            -c {params.clust[coverage]} \
            --cov-mode {params.clust[cov_mode]}
        rm -rf "$CHUNK_TMP"

        mmseqs createsubdb "$CHUNK_CLUSTERDB" "$CHUNK_DB" "$REPDB"
        mmseqs createsubdb "$CHUNK_CLUSTERDB" "${{CHUNK_DB}}_h" "${{REPDB}}_h"
        mmseqs createtsv "$CHUNK_DB" "$CHUNK_DB" "$CHUNK_CLUSTERDB" "{output.clusters}"

        verify_mmseqs_db "$REPDB"
        verify_mmseqs_db "${{REPDB}}_h"

        rm -rf "$BUILD_DIR"
        """

rule centroid_merge_querydb:
    input:
        chunk_repdbs = expand(
            "data/01.sOTUs_centroid/{mode}/chunks/{chunk}.repDB.dbtype",
            chunk=CENTROID_CHUNK_IDS, allow_missing=True
        ),
        chunk_clusters = expand(
            "data/01.sOTUs_centroid/{mode}/chunks/{chunk}.clusters.tsv",
            chunk=CENTROID_CHUNK_IDS, allow_missing=True
        ),
        exact_dups = "data/01.sOTUs_centroid/{mode}/dedup/exact_dups.tsv"
    output:
        db = "data/01.sOTUs_centroid/{mode}/queryDB.dbtype",
        clusters = "data/01.sOTUs_centroid/{mode}/queryDB_clusters.tsv"
    threads: 16
    resources:
        local_mem_mb=128000,
        runtime=2880
    conda: "../" + config["tools"]["DISCOVERY_ENV"]
    params:
        clust = lambda wildcards: config["querydb_clustering"],
        chunk_ids = CENTROID_CHUNK_IDS
    shell:
        """
        set -euo pipefail
        export LC_ALL=C

        verify_mmseqs_db() {{
            local db="$1" tries=0
            while [ ! -s "$db.dbtype" ] || [ ! -s "$db.index" ]; do
                tries=$((tries + 1))
                if [ "$tries" -ge 10 ]; then
                    echo "ERROR: $db is missing or incomplete after ${{tries}} checks over ~20s." >&2
                    return 1
                fi
                sleep 2
            done
        }}

        MODE_DIR="data/01.sOTUs_centroid/{wildcards.mode}"
        CHUNK_DIR="$MODE_DIR/chunks"
        BUILD_DIR="{resources.tmpdir}/centroid_merge_{wildcards.mode}"
        rm -rf "$BUILD_DIR" "$MODE_DIR"/queryDB*
        mkdir -p "$BUILD_DIR"

        AVAIL_MEM_MB=$(free -m | awk '/^Mem:/{{print $7}}')
        if [ -z "$AVAIL_MEM_MB" ] || [ "$AVAIL_MEM_MB" -le 0 ]; then
            AVAIL_MEM_MB={resources.local_mem_mb}
        fi
        SORT_BUF_MB=$(( AVAIL_MEM_MB / 2 ))
        SORT_BUF_CAP_MB=64000
        if [ "$SORT_BUF_MB" -gt "$SORT_BUF_CAP_MB" ]; then
            SORT_BUF_MB=$SORT_BUF_CAP_MB
        fi

        cat {input.chunk_clusters} > "$BUILD_DIR/chunk_clusters.tsv"

        # join-based expand (not awk-array) for why an in-memory array OOMs
        # here.
        sort -t $'\t' -k1,1 -S "${{SORT_BUF_MB}}M" -T "$BUILD_DIR" \
            --parallel={threads} --compress-program=pigz \
            "{input.exact_dups}" > "$BUILD_DIR/exact_dups.by_survivor.tsv"
        sort -t $'\t' -k2,2 -S "${{SORT_BUF_MB}}M" -T "$BUILD_DIR" \
            --parallel={threads} --compress-program=pigz \
            "$BUILD_DIR/chunk_clusters.tsv" > "$BUILD_DIR/chunk_clusters.by_member.tsv"

        {{
            cat "$BUILD_DIR/chunk_clusters.tsv"
            join -t $'\t' -1 1 -2 2 -o 2.1,1.2 \
                "$BUILD_DIR/exact_dups.by_survivor.tsv" "$BUILD_DIR/chunk_clusters.by_member.tsv"
        }} > "$BUILD_DIR/chunk_expanded.tsv"

        rm -f "$BUILD_DIR/exact_dups.by_survivor.tsv" "$BUILD_DIR/chunk_clusters.by_member.tsv"

        # Balanced pairwise merge tree over the 16 chunk repDBs; concatdbs
        # only takes two DBs at a time and must stay --threads 1, so
        # independent pairs within a round run backgrounded.
        CURRENT=()
        for c in {params.chunk_ids}; do
            CURRENT+=("$CHUNK_DIR/$c.repDB")
        done
        ROUND=0
        while [ "${{#CURRENT[@]}}" -gt 1 ]; do
            NEXT_ROUND=()
            PIDS=()
            for ((i = 0; i < ${{#CURRENT[@]}}; i += 2)); do
                if [ $((i + 1)) -lt "${{#CURRENT[@]}}" ]; then
                    A="${{CURRENT[$i]}}"
                    B="${{CURRENT[$((i + 1))]}}"
                    MERGED="$BUILD_DIR/merge_r${{ROUND}}_$((i / 2)).db"
                    (
                        verify_mmseqs_db "$A"
                        verify_mmseqs_db "${{A}}_h"
                        verify_mmseqs_db "$B"
                        verify_mmseqs_db "${{B}}_h"
                        mmseqs concatdbs "$A" "$B" "$MERGED" --threads 1 -v 1
                        mmseqs concatdbs "${{A}}_h" "${{B}}_h" "${{MERGED}}_h" --threads 1 -v 1
                        verify_mmseqs_db "$MERGED"
                        verify_mmseqs_db "${{MERGED}}_h"
                    ) &
                    PIDS+=($!)
                    NEXT_ROUND+=("$MERGED")
                else
                    NEXT_ROUND+=("${{CURRENT[$i]}}")
                fi
            done
            for PID in "${{PIDS[@]}}"; do
                wait "$PID"
            done
            CURRENT=("${{NEXT_ROUND[@]}}")
            ROUND=$((ROUND + 1))
        done
        ACC="${{CURRENT[0]}}"
        rm -rf "$CHUNK_DIR"

        mv "$ACC" "$BUILD_DIR/combinedDB"
        mv "$ACC.dbtype" "$BUILD_DIR/combinedDB.dbtype"
        mv "$ACC.index" "$BUILD_DIR/combinedDB.index"
        [ -f "$ACC.source" ] && mv "$ACC.source" "$BUILD_DIR/combinedDB.source"
        mv "${{ACC}}_h" "$BUILD_DIR/combinedDB_h"
        mv "${{ACC}}_h.dbtype" "$BUILD_DIR/combinedDB_h.dbtype"
        mv "${{ACC}}_h.index" "$BUILD_DIR/combinedDB_h.index"

        mmseqs linclust "$BUILD_DIR/combinedDB" "$BUILD_DIR/clusterDB" "$BUILD_DIR/clust_tmp" \
            --threads {threads} \
            --min-seq-id {params.clust[min_seq_id]} \
            -c {params.clust[coverage]} \
            --cov-mode {params.clust[cov_mode]}

        mmseqs createsubdb "$BUILD_DIR/clusterDB" "$BUILD_DIR/combinedDB" "$MODE_DIR/queryDB"
        mmseqs createsubdb "$BUILD_DIR/clusterDB" "$BUILD_DIR/combinedDB_h" "$MODE_DIR/queryDB_h"

        mmseqs createtsv "$BUILD_DIR/combinedDB" "$BUILD_DIR/combinedDB" "$BUILD_DIR/clusterDB" "$BUILD_DIR/combined_clusters.tsv"

        sort -t $'\t' -k2,2 -S "${{SORT_BUF_MB}}M" -T "$BUILD_DIR" \
            --parallel={threads} --compress-program=pigz \
            "$BUILD_DIR/combined_clusters.tsv" > "$BUILD_DIR/combined_clusters.sorted.tsv"
        sort -t $'\t' -k1,1 -S "${{SORT_BUF_MB}}M" -T "$BUILD_DIR" \
            --parallel={threads} --compress-program=pigz \
            "$BUILD_DIR/chunk_expanded.tsv" > "$BUILD_DIR/chunk_expanded.sorted.tsv"

        join -t $'\t' -1 2 -2 1 -o 1.1,2.2 \
            "$BUILD_DIR/combined_clusters.sorted.tsv" "$BUILD_DIR/chunk_expanded.sorted.tsv" \
            > "{output.clusters}"

        rm -rf "$BUILD_DIR"
        """

rule mmseqs_build_refdb:
    input: "data/01.sOTUs_centroid/{mode}/ncbi_nucleotide_reference.fasta"
    output: "data/01.sOTUs_centroid/{mode}/refDB.dbtype"
    run:
        prefix = os.path.join("data/01.sOTUs_centroid/{mode}".format(mode=wildcards.mode), "refDB")
        shell(f"rm -rf {prefix}* && mmseqs createdb {input} {prefix}")

rule mmseqs_search:
    input:
        query = "data/01.sOTUs_centroid/{mode}/queryDB.dbtype",
        ref   = "data/01.sOTUs_centroid/{mode}/refDB.dbtype"
    output:
        hits = "data/01.sOTUs_centroid/{mode}/hits_raw.tsv"
    params:
        q_prefix = "data/01.sOTUs_centroid/{mode}/queryDB",
        r_prefix = "data/01.sOTUs_centroid/{mode}/refDB",
        mmseqs   = lambda wildcards: config["scenarios"][wildcards.mode]["mmseqs_params"]
    threads: 16
    # No mem_mb: this cluster's `shared` partition schedules on cores
    # only and rejects any --mem request outright.
    resources:
        local_mem_mb=128000,
        runtime=2880
    shell:
        """
        TMP_DIR="{resources.tmpdir}/mmseqs_search_{wildcards.mode}"
        rm -rf "$TMP_DIR" && mkdir -p "$TMP_DIR"
        sync; sleep 5

        mmseqs search "{params.q_prefix}" "{params.r_prefix}" "$TMP_DIR/res" "$TMP_DIR" \
            --threads {threads} \
            --search-type 2 \
            --translate 0 \
            --split-mode 2 \
            -c "{params.mmseqs[coverage_threshold]}" \
            --min-seq-id "{params.mmseqs[min_seq_id]}"

        sync

        mmseqs convertalis "{params.q_prefix}" "{params.r_prefix}" "$TMP_DIR/res" "{output.hits}" \
            --threads {threads} --format-output query,target,pident,alnlen,evalue,fident

        rm -rf "$TMP_DIR"
        """

# Best-hit max_seq_id/min_seq_id enforcement step 15): convertalis can report
# several partial alignments per query, so "best hit" means the single top-
# scoring alignment (lowest evalue), not the maximum fident seen anywhere for
# that query -- a short, spuriously high-identity fragment must not win over a
# longer, better-scoring real alignment. fident is dropped from the output so
# expand_cluster_hits/build_matrix/extract_sra_ids need no changes.
rule filter_best_hit:
    input:
        hits_raw = "data/01.sOTUs_centroid/{mode}/hits_raw.tsv"
    output:
        "data/01.sOTUs_centroid/{mode}/hits_filtered.tsv"
    params:
        mmseqs = lambda wildcards: config["scenarios"][wildcards.mode]["mmseqs_params"]
    shell:
        """
        if [ ! -s "{input.hits_raw}" ]; then
            : > "{output[0]}"
            exit 0
        fi
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/filter_best_hit.py \
            --hits-raw "{input.hits_raw}" \
            --min-seq-id "{params.mmseqs[min_seq_id]}" \
            --max-seq-id "{params.mmseqs[max_seq_id]}" \
            --out "{output[0]}"
        """

# Credit every collapsed run behind each hit, not just whichever
# fragment mmseqs kept as cluster representative.
rule expand_cluster_hits:
    input:
        hits_raw = "data/01.sOTUs_centroid/{mode}/hits_filtered.tsv",
        clusters = "data/01.sOTUs_centroid/{mode}/queryDB_clusters.tsv"
    output:
        "data/01.sOTUs_centroid/{mode}/hits.tsv"
    threads: 8
    resources:
        local_mem_mb=64000,
        runtime=360,
    retries: 2
    shell:
        """
        set -euo pipefail
        export LC_ALL=C

        if [ ! -s "{input.hits_raw}" ]; then
            : > "{output[0]}"
            exit 0
        fi

        BUILD_DIR="{resources.tmpdir}/expand_cluster_hits_{wildcards.mode}"
        rm -rf "$BUILD_DIR"
        mkdir -p "$BUILD_DIR"

        AVAIL_MEM_MB=$(free -m | awk '/^Mem:/{{print $7}}')
        if [ -z "$AVAIL_MEM_MB" ] || [ "$AVAIL_MEM_MB" -le 0 ]; then
            AVAIL_MEM_MB={resources.local_mem_mb}
        fi
        SORT_BUF_MB=$(( AVAIL_MEM_MB / 2 ))
        SORT_BUF_CAP_MB=32000
        if [ "$SORT_BUF_MB" -gt "$SORT_BUF_CAP_MB" ]; then
            SORT_BUF_MB=$SORT_BUF_CAP_MB
        fi

        # queryDB_clusters.tsv can be huge (one row per raw fragment), so
        # only pull membership rows for reps that hit something.
        cut -f1 "{input.hits_raw}" | sort -u \
            -S "${{SORT_BUF_MB}}M" -T "$BUILD_DIR" --parallel={threads} \
            > "$BUILD_DIR/reps.tsv"
        awk -F'\t' 'NR==FNR {{keep[$1]=1; next}} ($1 in keep)' \
            "$BUILD_DIR/reps.tsv" "{input.clusters}" > "$BUILD_DIR/members.tsv"

        # Same sorted-merge `join` idiom as centroid_merge_querydb's
        # first/second-expand fixes: an in-memory array/dataframe keyed by
        # cluster rep isn't bounded by any sort buffer and OOMs regardless of
        # available RAM, since this cluster's shared nodes have no SLURM-level
        # memory enforcement.
        sort -t $'\t' -k1,1 -S "${{SORT_BUF_MB}}M" -T "$BUILD_DIR" \
            --parallel={threads} --compress-program=pigz \
            "{input.hits_raw}" > "$BUILD_DIR/hits.sorted.tsv"
        sort -t $'\t' -k1,1 -S "${{SORT_BUF_MB}}M" -T "$BUILD_DIR" \
            --parallel={threads} --compress-program=pigz \
            "$BUILD_DIR/members.tsv" > "$BUILD_DIR/members.sorted.tsv"

        # Left outer join on query/rep id; no membership row for a rep
        # (shouldn't normally happen) falls back to crediting the rep
        # itself rather than dropping the hit.
        join -t $'\t' -1 1 -2 1 -a 1 -e NOMATCH \
            -o 1.1,1.2,1.3,1.4,1.5,2.2 \
            "$BUILD_DIR/hits.sorted.tsv" "$BUILD_DIR/members.sorted.tsv" \
            > "$BUILD_DIR/joined.tsv"

        awk -F'\t' 'BEGIN{{OFS="\t"}} {{member=($6=="NOMATCH")?$1:$6; print member,$2,$3,$4,$5}}' \
            "$BUILD_DIR/joined.tsv" \
            | sort -u -t $'\t' -k1,1 -k2,2 \
                -S "${{SORT_BUF_MB}}M" -T "$BUILD_DIR" --parallel={threads} \
            > "{output[0]}"

        rm -rf "$BUILD_DIR"
        """

rule build_matrix:
    input: hits="data/01.sOTUs_centroid/{mode}/hits.tsv", unanchored="data/01.sOTUs_centroid/{mode}/unanchored_palm_ids.txt"
    output: "data/01.sOTUs_centroid/{mode}/unique.csv"
    params:
        db = config["database"]
    script:
        "../scripts/align_matrix.py"

rule extract_sra_ids:
    input:
        hits="data/01.sOTUs_centroid/{mode}/hits.tsv"
    output:
        "data/01.sOTUs_centroid/{mode}/sra_ids.txt"
    script:
        "../scripts/extract_sra_ids.py"
