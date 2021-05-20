rule samtools_sort:
    input:
        expand("results/{library}/{library}_{type}.sam", type=["mutated", "base"], allow_missing=True)
    output:
        temp(expand("results/{library}/{library}_{type}.bam", type=["mutated", "base"], allow_missing=True))
    log:
        "results/logs/samtools_sort_{library}.log"
    benchmark:
        "results/benchmarks/samtools_sort_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools sort --threads {threads} {input[0]} > {output[0]} ; "
        "samtools sort --threads {threads} {input[1]} > {output[1]}"
