rule circos_cov:
    input:
        "results/reference/{reference}.fa.fai",
        expand("results/alns/{library}_on_{reference}.regions.bed", library=config["samples"], allow_missing=True)
    output:
        "results/alns/circos_cov_{reference}.png"
    log:
        "results/logs/circos_cov_{reference}.log"
    benchmark:
        "results/benchmarks/circos_cov_{reference}.benchmark.txt"
    params:
        reference="{reference}",
        library=config["samples"]
    threads: 1
    resources:
        mem_mb=1000
    script:
        "../scripts/circos_cov.R"
