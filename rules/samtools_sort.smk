sample=[config["base"], config["mutated"]]

rule samtools_sort:
    input:
        "results/alignments/bwa/{sample}.sam"
    output:
        "results/alignments/bwa/{sample}.bam"
    log:
        "results/logs/samtools_sort_{sample}.log"
    benchmark:
        "results/benchmarks/samtools_sort_{sample}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools sort --threads {threads} {input} > {output}"

# samtools sort {input} > {output}
# samtools sort --threads {threads} {input.bam} -m 2G - > ${idSample}_${idRun}.bam