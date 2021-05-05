sample=expand("{file}{type}", file=[config["base"], config["mutated"]], type=["", ".md"]) # can add BQSR and all needed here

rule samtools_index:
    input:
        "results/alignments/bwa/{sample}.bam"
    output:
        "results/alignments/bwa/{sample}.bam.bai"
    log:
        "results/logs/samtools_index_{sample}.log"
    benchmark:
        "results/benchmarks/samtools_index_{sample}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input}"
