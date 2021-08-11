rule strelka2sql:
    input:
        expand("results/mutations/{vcfs}_strelka2.tsv", vcfs=config["vcfs"])
    output:
        temp("results/strelka2_raw.csv"),
        "results/strelka2_raw.sql"
    log:
        "results/logs/strelka2sql.log"
    benchmark:
        "results/benchmarks/strelka2sql.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/strelka2sql.R"
        