rule gatk_haplotypecaller:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        "results/{library}/{library}_mutated.md.bam",
        "results/{library}/{library}_mutated.md.bam.bai",
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        "results/{library}/gatk/{library}.g.vcf"
    log:
        "results/logs/gatk_haplotypecaller_{library}.log"
    benchmark:
        "results/benchmarks/gatk_haplotypecaller_{library}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk HaplotypeCaller -R {input[0]} -I {input[1]} -O {output} -ERC GVCF"