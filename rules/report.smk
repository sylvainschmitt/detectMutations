rule report:
    input:
        index=expand("results/reference/{reference}.fa.fai", reference=config["references"]),
        coverage=expand("results/alns/{library}_on_{reference}.regions.bed", library=config["samples"], reference=config["references"]),
        circos_cov=expand("results/alns/circos_cov_{reference}.png", reference=config["references"]),
        mutations="results/mutations/mutations.tsv",
        cv="results/mutations/cross_validation.tsv",
        spectra="results/mutations/spectra.tsv"
    output:
        "results/report.html"
    log:
        "results/logs/report.log"
    benchmark:
        "results/benchmarks/report.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    script:
        "../report/report.Rmd"
