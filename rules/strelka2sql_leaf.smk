rule strelka2sql_leaf:
    input:
        expand("results/mutations_leaf/{tumor}_vs_{base}.nonhz.tsv", tumor=config["leaf"], base=config["cambium_ref"])
    output:
        temp("results/mutations_leaf/leaf_nonhz_mutations.csv"),
        "results/mutations_leaf/leaf_nonhz_mutations.sql"
    log:
        "results/logs/strelka2sql_leaf.log"
    benchmark:
        "results/benchmarks/strelka2sql_leaf.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/strelka2sql.R"
