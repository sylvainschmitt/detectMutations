rule psl2pos:
    input:
        "results/napoleon/napoleon_mutations.tsv",
        "results/napoleon/napoleon_mutations.psl"
    output:
        "results/napoleon_mutations.tsv"
    log:
        "results/logs/psl2pos.log"
    benchmark:
        "results/benchmarks/psl2pos.benchmark.txt"
    script:
        "../scripts/psl2pos.R"
