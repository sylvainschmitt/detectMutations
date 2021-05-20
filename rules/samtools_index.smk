rule samtools_index:
    input:
        expand("results/{library}/{library}_{type}.bam", type=["mutated", "base"], allow_missing=True)
    output:
        temp(expand("results/{library}/{library}_{type}.bam.bai", type=["mutated", "base"], allow_missing=True))
    log:
        "results/logs/samtools_index_{library}.log"
    benchmark:
        "results/benchmarks/samtools_index_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input[0]} ; samtools index {input[1]}"
