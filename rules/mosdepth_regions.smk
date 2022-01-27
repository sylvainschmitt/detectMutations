rule mosdepth_regions:
    input:
        "results/reference/{reference}.fa",
        "results/alns/{library}_on_{reference}.md.cram",
        expand("results/reference/{reference}.{ext}", 
                ext=["fa.amb", "fa.ann", "fa.bwt", "fa.pac", "fa.sa"], allow_missing=True),
        "results/alns/{library}_on_{reference}.md.cram.crai"
    output:
        "results/alns/{library}_on_{reference}.mosdepth.global.dist.txt",
        "results/alns/{library}_on_{reference}.mosdepth.region.dist.txt"
    log:
        "results/logs/mosdepth_regions_{library}_{reference}.log"
    benchmark:
        "results/benchmarks/mosdepth_regions_{library}_{reference}.benchmark.txt"
    singularity: 
        "docker://quay.io/biocontainers/mosdepth:0.2.4--he527e40_0"
    shell:
        "mosdepth -n --fast-mode --by 10000 -t {threads} -f {input[0]} results/alns/{wildcards.library} {input[1]}"
