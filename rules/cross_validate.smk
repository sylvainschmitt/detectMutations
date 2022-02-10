rule cross_validate:
    input:
        expand("results/mutations/{reference}_mutations_on_{references}.tsv", reference=config["reference"], references=config["references"])
    output:
        "results/mutations/cross_validation.tsv"
    log:
        "results/logs/cross_validate.log"
    benchmark:
        "results/benchmarks/cross_validate.benchmark.txt"
    script:
        "../scripts/cross_validate.R"
