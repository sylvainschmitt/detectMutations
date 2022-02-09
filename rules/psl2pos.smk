rule psl2pos:
    input:
        "results/mutations/mutations.tsv",
        expand("results/mutations/{reference}_mutations_on_{references}.psl", reference=config["reference"], allow_missing=True)
    output:
        expand("results/mutations/{reference}_mutations_on_{references}.tsv", reference=config["reference"], allow_missing=True)
    log:
        "results/logs/psl2pos_{references}.log"
    benchmark:
        "results/benchmarks/psl2pos.benchmark_{references}.txt"
    script:
        "../scripts/psl2pos.R"
