rule tsv2sql:
    input:
        expand("results/mutations/{vcfs}_on_{chromosome}_{caller}.tsv", vcfs=config["vcfs"], chromosome=chromosomes, allow_missing=True)
    output:
        temp("results/{caller}_raw.csv"),
        "results/{caller}_raw.sql"
    log:
        "results/logs/{caller}_tsv2sql.log"
    benchmark:
        "results/benchmarks/{caller}_tsv2sql.benchmark.txt"
    script:
        "../scripts/tsv2sql.R"
