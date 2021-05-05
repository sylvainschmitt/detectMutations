sample=[config["base"]]

rule gatk_genotypegvcfs:
    input:
        "results/germline/gatk/{sample}.g.vcf",
         multiext("results/raw_data/reference/reference", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict")
    output:
        "results/germline/gatk/{sample}.vcf"
    log:
        "results/logs/gatk_genotypegvcfs_{sample}.log"
    benchmark:
        "results/benchmarks/gatk_genotypegvcfs_{sample}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk GenotypeGVCFs -R {input[1]} -V {input[0]} -O {output}"
        