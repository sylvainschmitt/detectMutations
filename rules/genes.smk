rule genes:
    input:
        "results/mutations/mutations_filtered.tsv",
        expand("{refdir}{genes}", refdir=config["refdir"], genes=config["genes"]),
        expand("{refdir}{annotation}", refdir=config["refdir"], annotation=config["annotation"])
    output:
        "results/mutations/genes.tsv"
    log:
        "results/logs/genes.log"
    benchmark:
        "results/benchmarks/genes.benchmark.txt"
    script:
        "../scripts/genes.R"
