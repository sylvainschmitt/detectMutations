rule bcftools_vcf:
    input:
        "results/variants/{reference}.mpileup"
    output:
        "results/variants/{reference}.vcf"
    log:
        "results/logs/bcftools_vcf_{reference}.log"
    benchmark:
        "results/benchmarks/bcftools_vcf_{reference}.benchmark.txt"
    shell:
        "bcftools call -v -m {input[0]} > {output[0]}"
