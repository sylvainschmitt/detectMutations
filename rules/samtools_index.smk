rule samtools_index:
    input:
        "results/{library}/{library}.sorted.cram"
    output:
        temp("results/{library}/{library}.sorted.cram.crai")
    log:
        "results/logs/samtools_index_{library}.log"
    benchmark:
        "results/benchmarks/samtools_index_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input}"
