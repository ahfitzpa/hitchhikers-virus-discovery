# One-time database provisioning, driven by setup_databases.sh (see SETUP.md).
# Every rule is idempotent -- it checks for its own.done flag or target files
# before downloading anything. The MMseqs2 nt and ICCTax setup rules are gone
# three-tool taxonomy switch.

rule setup_kraken2:
    output:
        hash=os.path.join(config["databases"]["kraken2_db"], "hash.k2d"),
        ktax=os.path.join(config["databases"]["kraken2_db"], "ktaxonomy.tsv")
    params:
        url=config["kraken2_db_url"]
    benchmark:
        "benchmarks/setup_kraken2/setup.tsv"
    shell:
        """
        DB_DIR=$(dirname {output.hash})
        mkdir -p "$DB_DIR"
        wget -q -O "$DB_DIR/k2_db.tar.gz" {params.url}
        tar -xzf "$DB_DIR/k2_db.tar.gz" -C "$DB_DIR"
        rm -f "$DB_DIR/k2_db.tar.gz"
        """

rule setup_genomad:
    conda: get_env("GENOMAD_ENV")
    output:
        version=os.path.join(config["databases"]["genomad_db"], "version.txt")
    benchmark:
        "benchmarks/setup_genomad/setup.tsv"
    shell:
        # download-database writes a genomad_db/ subdir under the path it's
        # given, so point it at the configured path's parent.
        """
        PARENT=$(dirname {config[databases][genomad_db]})
        mkdir -p "$PARENT"
        genomad download-database "$PARENT"
        """

rule setup_virsorter2:
    conda: get_env("TAX_QC_ENV")
    output:
        dir=directory(config["databases"]["virsorter2_db"]),
        flag=os.path.join(config["databases"]["virsorter2_db"], ".done")
    benchmark:
        "benchmarks/setup_virsorter2/setup.tsv"
    shell:
        """
        if [ ! -f "{output.flag}" ]; then
            virsorter setup -d {output.dir} -j 4
            virsorter config --init-source --db-dir={output.dir}
        fi
        touch {output.flag}
        """

rule setup_checkv:
    conda: get_env("TAX_QC_ENV")
    output:
        dir=directory(config["databases"]["checkv_db"]),
        flag=os.path.join(config["databases"]["checkv_db"], ".done")
    benchmark:
        "benchmarks/setup_checkv/setup.tsv"
    shell:
        """
        if [ ! -f "{output.flag}" ]; then
            checkv download_database {output.dir}
        fi
        touch {output.flag}
        """

rule setup_rdrpcatch:
    conda: get_env("TAX_CORE_ENV")
    output:
        dir=directory(config["databases"]["rdrp_catch_hmm"]),
        flag=os.path.join(config["databases"]["rdrp_catch_hmm"], ".done")
    benchmark:
        "benchmarks/setup_rdrpcatch/setup.tsv"
    shell:
        """
        if [ ! -f "{output.flag}" ]; then
            rdrpcatch databases --destination-dir {output.dir}
        fi
        touch {output.flag}
        """

rule setup_metacerberus:
    output:
        sentinel = touch(os.path.join(config["databases"]["metacerberus"], ".done"))
    conda: get_env("META_ENV")
    params:
        dbs = " ".join(config["annotation"]["metacerberus_download"])
    benchmark:
        "benchmarks/setup_metacerberus/setup.tsv"
    shell:
        # Only the HMM sets run_metacerberus actually uses, not --download's
        # all-databases default.
        """
        mkdir -p {config[databases][metacerberus]}
        metacerberus.py --setup
        metacerberus.py --download {params.dbs} --db-path {config[databases][metacerberus]}
        """
