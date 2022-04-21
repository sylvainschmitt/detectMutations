rule join_all:
    input:
        mutations="results/mutations/mutations_filtered.tsv",
        validation="results/mutations/cross_validation.tsv",
        spectra="results/mutations/spectra.tsv",
        genes="results/mutations/genes.tsv",
        te="results/mutations/te.tsv"
    output:
        "results/mutations.tsv"
    log:
        "results/logs/join_all.log"
    benchmark:
        "results/benchmarks/join_all.benchmark.txt"
    script:
        "../scripts/join_all.R"
