rule gatk_genotypegvcfs:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        "results/{library}/gatk/{library}.g.vcf",
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        "results/{library}/gatk/{library}.unfiltered.vcf"
    log:
        "results/logs/gatk_genotypegvcfs_{library}.log"
    benchmark:
        "results/benchmarks/gatk_genotypegvcfs_{library}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk GenotypeGVCFs -R {input[0]} -V {input[1]} -O {output}"
        