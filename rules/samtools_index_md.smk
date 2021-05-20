rule samtools_index_md:
    input:
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True)
    output:
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True)
    log:
        "results/logs/samtools_index_md_{library}.log"
    benchmark:
        "results/benchmarks/samtools_index_md_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input[0]} ; samtools index {input[1]}"
