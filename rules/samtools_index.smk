rule samtools_index:
    input:
        "results/{library}/{library}.bam"
    output:
        temp("results/{library}/{library}.bam.bai")
    log:
        "results/logs/samtools_index_{library}.log"
    benchmark:
        "results/benchmarks/samtools_index_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input}"
