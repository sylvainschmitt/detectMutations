rule bedtools_makewindows:
    input:
        expand("results/reference/{reference}.fa.fai", reference=config["reference"])
    output:
        temp(expand("results/reference/{reference}.1kb.bed", reference=config["reference"]))
    log:
        "results/logs/bedtools_makewindows.log"
    benchmark:
        "results/benchmarks/bedtools_makewindows.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "cut -f 1,2 {input} | bedtools makewindows -g stdin -w 10000 > {output}"
        