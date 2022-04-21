rule te:
    input:
        "results/mutations/mutations_filtered.tsv",
        expand("{refdir}{te}", refdir=config["refdir"], te=config["te"])
    output:
        "results/mutations/te.tsv"
    log:
        "results/logs/te.log"
    benchmark:
        "results/benchmarks/te.benchmark.txt"
    script:
        "../scripts/te.R"
