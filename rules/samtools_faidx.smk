rule samtools_faidx:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"])
    output:
        expand("results/reference/{reference}.fa.fai", reference=config["reference"])
    log:
        "results/logs/samtools_faidx.log"
    benchmark:
        "results/benchmarks/samtools_faidx.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools faidx {input}"
