rule spectra:
    input:
        "results/mutations/mutations_filtered.tsv",
         expand("results/reference/{reference}.fa", reference=config["reference"]),
         expand("results/reference/{reference}.fa.fai", reference=config["reference"])
    output:
        "results/mutations/spectra.tsv"
    log:
        "results/logs/spectra.log"
    benchmark:
        "results/benchmarks/spectra.benchmark.txt"
    script:
        "../scripts/spectra.R"
