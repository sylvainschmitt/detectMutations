rule filter_mutations:
    input:
        "results/mutations/mutations.raw.sql"
    output:
        "results/mutations/mutations.tsv",
    log:
        "results/logs/filter_mutations.log"
    benchmark:
        "results/benchmarks/filter_mutations.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    script:
        "../scripts/filter_mutations.R"
