rule samtools_faidx:
    input:
        expand("results/reference/{reference}_{chromosome}.fa", reference=config["reference"], allow_missing = True)
    output:
        expand("results/reference/{reference}_{chromosome}.fa.fai", reference=config["reference"], allow_missing = True)
    log:
        "results/logs/samtools_faidx_{chromosome}.log"
    benchmark:
        "results/benchmarks/samtools_faidx_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools faidx {input}"
