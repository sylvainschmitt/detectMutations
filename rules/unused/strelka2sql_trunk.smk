rule strelka2sql_trunk:
    input:
        expand("results/mutations/{trunk}_vs_{base}.raw.tsv", trunk=config["trunk"], base=config["base"])
    output:
        temp("results/trunk_raw_mutations.csv"),
        "results/trunk_raw_mutations.sql"
    log:
        "results/logs/strelka2sql_trunk.log"
    benchmark:
        "results/benchmarks/strelka2sql_trunk.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/strelka2sql.R"
