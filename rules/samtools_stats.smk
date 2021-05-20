rule samtools_stats:
    input:
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True)
    output:
        temp(expand("results/{library}/{library}_{type}.md.bam.stats.out", type=["mutated", "base"], allow_missing=True))
    log:
        "results/logs/samtools_stats_{library}.log"
    benchmark:
        "results/benchmarks/samtools_stats_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools stats --threads {threads} {input[0]} > {output[0]} ; "
        "samtools stats --threads {threads} {input[1]} > {output[1]}"
