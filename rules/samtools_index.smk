rule samtools_index:
    input:
        "results/{library}/{library}_{chromosome}.sorted.cram"
    output:
        temp("results/{library}/{library}_{chromosome}.sorted.cram.crai")
    log:
        "results/logs/samtools_index_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/samtools_index_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input}"
