rule napomutations2bed:
    input:
        "results/napoleon/napoleon_mutations.tsv",
        "results/napoleon/napoleon.fa.fai"
    output:
        temp("results/napoleon/napoleon_mutations.bed")
    log:
        "results/logs/napomutations2bed.log"
    benchmark:
        "results/benchmarks/napomutations2bed.benchmark.txt"
    script:
        "../scripts/napomutations2bed.R"
