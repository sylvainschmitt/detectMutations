rule samtools_faidx:
    input:
        "results/raw_data/reference/reference.fa"
    output:
        "results/raw_data/reference/reference.fa.fai"
    log:
        "results/logs/samtools_faidx.log"
    benchmark:
        "results/benchmarks/samtools_faidx.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools faidx {input}"
