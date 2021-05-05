sample=[config["base"], config["mutated"]]

rule samtools_view:
    input:
        "results/alignments/bwa/{sample}.bam"
    output:
        "results/alignments/bwa/{sample}.paired.bam"
    log:
        "results/logs/samtools_view_{sample}.log"
    benchmark:
        "results/benchmarks/samtools_view_{sample}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools view -u -f 1 -F 12 {input} > {output}"
