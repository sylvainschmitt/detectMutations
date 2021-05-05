sample=[config["base"]]

rule gatk_haplotypecaller:
    input:
        "results/alignments/bwa/{sample}.md.bam",
         multiext("results/raw_data/reference/reference", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict")
    output:
        "results/germline/gatk/{sample}.g.vcf"
    log:
        "results/logs/gatk_haplotypecaller_{sample}.log"
    benchmark:
        "results/benchmarks/gatk_haplotypecaller_{sample}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk HaplotypeCaller -R {input[1]} -I {input[0]} -O {output} -ERC GVCF"
