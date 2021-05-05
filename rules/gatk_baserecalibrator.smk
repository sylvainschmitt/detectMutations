sample=[config["base"], config["mutated"]]

rule gatk_baserecalibrator:
    input:
        "results/alignments/bwa/{sample}.md.bam",
         multiext("results/raw_data/reference/reference", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict")
    output:
        "results/alignments/bwa/{sample}.recal.table"
    log:
        "results/logs/gatk_baserecalibrator_{sample}.log"
    benchmark:
        "results/benchmarks/gatk_baserecalibrator_{sample}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk BaseRecalibrator -I {input[0]} -R {input[1]} -O {output}"
        