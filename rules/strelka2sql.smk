rule strelka2sql:
    input:
        "results/mutations/{tumor}_vs_{base}.raw.tsv"
    output:
        temp("results/mutations/{tumor}_vs_{base}.raw.csv"),
        "results/mutations/{tumor}_vs_{base}.raw.sql"
    log:
        "results/logs/strelka2sql_{tumor}_vs_{base}.log"
    benchmark:
        "results/benchmarks/strelka2sql_{tumor}_vs_{base}.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/strelka2sql.R"
