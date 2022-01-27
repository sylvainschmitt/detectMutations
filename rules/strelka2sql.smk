rule strelka2sql:
    input:
        expand("results/mutations/{comp}_on_{reference}.raw.tsv", comp=config["comps"], reference=config["references"])
    output:
        temp("results/mutations/{mutations.raw.csv"),
        "results/mutations/mutations.raw.sql"
    log:
        "results/logs/strelka2sql.log"
    benchmark:
        "results/benchmarks/strelka2sql.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/strelka2sql.R"
