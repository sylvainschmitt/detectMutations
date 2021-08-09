rule mosdepth:
    input:
        "results/{library}/{library}.md.cram",
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        "results/{library}/{library}.md.cram.crai"
    output:
        "results/{library}/{library}.regions.bed.gz",
        temp("results/{library}/{library}.mosdepth.global.dist.txt"),
        temp("results/{library}/{library}.mosdepth.region.dist.txt"),
        temp("results/{library}/{library}.regions.bed.gz.csi")
    log:
        "results/logs/mosdepth_{library}.log"
    benchmark:
        "results/benchmarks/mosdepth_{library}.benchmark.txt"
    singularity: 
        "docker://quay.io/biocontainers/mosdepth:0.2.4--he527e40_0"
    threads: 4
    shell:
        "mosdepth -t {threads} -n --fast-mode --by 10000 -f {input[1]} {wildcards.library} {input[0]}"
