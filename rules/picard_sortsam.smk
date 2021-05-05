sample=[config["base"], config["mutated"]]

rule picard_sortsam:
    input:
        "results/alignments/bwa/{sample}.sam"
    output:
        "results/alignments/bwa/{sample}.bam"
    log:
        "results/logs/picard_sortsam_{sample}.log"
    benchmark:
        "results/benchmarks/picard_sortsam_{sample}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/picard"
    shell:
        "java -Xmx4g -jar $PICARD SortSam I={input} O={output} SORT_ORDER=coordinate "
