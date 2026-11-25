rule fetch_sra_metadata:
    input:
        run_ids="data/01.sOTUs_centroid/{mode}/sra_ids.txt"
    output:
        filtered_csv="data/02.SRA_metadata/{mode}/filtered_unique.csv",
        metadata_tsv_gz="data/02.SRA_metadata/{mode}/unified_metadata.tsv.gz"
    conda: "../" + config["tools"]["DISCOVERY_ENV"]
    params:
        db_host = config["database"].get("host"),
        db_name = config["database"].get("name"),
        db_user = config["database"].get("user"),
        db_pass = config["database"].get("password")
    shell:
        """
        mkdir -p data/02.SRA_metadata/{wildcards.mode}
        
        export DB_HOST="{params.db_host}"
        export DB_NAME="{params.db_name}"
        export DB_USER="{params.db_user}"
        export PGPASSWORD="{params.db_pass}"
        
        # Use workflow.basedir to find the script reliably
        bash {workflow.basedir}/scripts/03.metadata_fetch.sh \
              "{wildcards.mode}" \
              "{input.run_ids}" \
              "{output.filtered_csv}" \
              "{output.metadata_tsv_gz}"
        """