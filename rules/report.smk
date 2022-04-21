rule report:
    input:
        scheme="dag/scheme.png",
        index=expand("results/reference/{reference}.fa.fai", reference=config["references"]),
        coverage=expand("results/alns/{library}_on_{reference}.regions.bed", library=config["samples"], reference=config["references"]),
        mutations_filterd="results/mutations/mutations_filtered.tsv",
        mutations="results/mutations.tsv",
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
        "../scripts/report.Rmd"
