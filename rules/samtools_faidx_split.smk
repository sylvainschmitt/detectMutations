rule samtools_faidx_split:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"])
    output:
        expand("results/reference/{reference}_{chromosome}.fa", reference=config["reference"], allow_missing = True)
    log:
        "results/logs/samtools_faidx_split_{chromosome}.log"
    benchmark:
        "results/benchmarks/samtools_faidx_split_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools faidx {input} {wildcards.chromosome} -o {output}"
