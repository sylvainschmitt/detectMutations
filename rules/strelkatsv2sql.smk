rule strelkatsv2sql:
    input:
        expand("results/mutations/{vcfs}_on_{chromosome}_strelka2.tsv", vcfs=config["vcfs"], chromosome=chromosomes)
    output:
        temp("results/strelka2_raw.csv"),
        "results/strelka2_raw.sql"
    log:
        "results/logs/strelkatsv2sql.log"
    benchmark:
        "results/benchmarks/strelkatsv2sql.benchmark.txt"
    script:
        "../scripts/strelkatsv2sql.R"
