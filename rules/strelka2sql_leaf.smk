rule strelka2sql_leaf:
    input:
        expand("results/mutations/{leaf}_vs_{base}.nontrunk.tsv", leaf=config["leaf"], base=config["base"])
    output:
        temp("results/leaf_nontrunk_mutations.csv"),
        "results/leaf_nontrunk_mutations.sql"
    log:
        "results/logs/strelka2sql_leaf.log"
    benchmark:
        "results/benchmarks/strelka2sql_leaf.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/strelka2sql.R"
