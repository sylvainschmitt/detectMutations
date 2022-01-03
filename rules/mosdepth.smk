rule mosdepth:
    input:
        "results/alns/{library}.md.cram"
    output:
        "results/alns/{library}.mosdepth.global.dist.txt",
        "results/alns/{library}.mosdepth.summary.txt"
    log:
        "results/logs/mosdepth_{library}.log"
    benchmark:
        "results/benchmarks/mosdepth_{library}.benchmark.txt"
    singularity: 
        "docker://quay.io/biocontainers/mosdepth:0.2.4--he527e40_0"
    shell:
        "mosdepth -n --fast-mode -t {threads} {wildcards.library} {input}"
