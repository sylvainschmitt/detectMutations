sample=[config["base"], config["mutated"]]

rule gatk_markduplicates:
    input:
        multiext("results/alignments/bwa/{sample}", ".bam", ".bam.bai")
    output:
        "results/alignments/bwa/{sample}.md.bam",
        "results/alignments/bwa/{sample}.bam.metrics"
    log:
        "results/logs/gatk_markduplicates_{sample}.log"
    benchmark:
        "results/benchmarks/gatk_markduplicates_{sample}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk MarkDuplicates I={input[0]} O={output[0]} M={output[1]}"
