rule bcftools_biallelic:
    input:
        "results/{reference}.vcf"
    output:
        temp("results/{reference}_biallelic.vcf")
    log:
        "results/logs/bcftools_biallelic_{reference}.log"
    benchmark:
        "results/benchmarks/bcftools_biallelic_{reference}.benchmark.txt"
    shell:
        "bcftools view --max-alleles 2 --threads {threads} {input} > {output}"
