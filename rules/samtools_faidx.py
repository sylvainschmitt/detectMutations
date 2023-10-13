rule samtools_faidx:
    input:
        "results/reference/{reference}.fa"
    output:
        "results/reference/{reference}.fa.fai"
    log:
        "results/logs/samtools_faidx_{reference}.log"
    benchmark:
        "results/benchmarks/samtools_faidx_{reference}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools faidx {input}"
