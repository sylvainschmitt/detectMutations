sample=[config["base"], config["mutated"]]

rule samtools_stats:
    input:
        "results/alignments/bwa/{sample}.md.bam"
    output:
        "results/alignments/bwa/{sample}.md.bam.stats.out"
    log:
        "results/logs/samtools_stats_{sample}.log"
    benchmark:
        "results/benchmarks/samtools_stats_{sample}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools stats --threads {threads} {input} > {output}"
