rule gatk2sql:
    input:
        expand("results/mutations/{chromosome}_gatk.tsv", chromosome=chromosomes)
    output:
        temp("results/gatk_raw.csv"),
        "results/gatk_raw.sql"
    log:
        "results/logs/gatk2sql.log"
    benchmark:
        "results/benchmarks/gatk2sql.benchmark.txt"
    script:
        "../scripts/gatk2sql.R"