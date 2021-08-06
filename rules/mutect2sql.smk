rule mutect2sql:
    input:
        expand("results/mutations/{vcfs}_mutect2.tsv", vcfs=config["vcfs"])
    output:
        temp("results/mutect2_raw.csv"),
        "results/mutect2_raw.sql"
    log:
        "results/logs/mutect2sql.log"
    benchmark:
        "results/benchmarks/mutect2sql.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.1/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/mutect2sql.R"
        