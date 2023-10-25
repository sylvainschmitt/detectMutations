rule gatk_snps:
    input:
        "results/{reference}_biallelic.vcf"
    output:
        temp("results/{reference}_intermediate.vcf"),
        "results/{reference}_filtered.vcf",
        temp("results/{reference}_biallelic.vcf.idx"),
        temp("results/{reference}_intermediate.vcf.idx"),
        temp("results/{reference}_filtered.vcf.idx")
    log:
        "results/logs/gatk_snps_{reference}.log"
    benchmark:
        "results/benchmarks/gatk_snps_{reference}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=2000
    shell:
        "gatk IndexFeatureFile -I {input} ;"
        "gatk VariantFiltration -V {input} -filter-expression \"QUAL < 30.0 || QD < 2.0 || FS > 60.0 || SOR > 3.0\" --filter-name \"FAIL\" "
        "-O {output[0]} ;"
        "gatk SelectVariants -V {output[0]} -select-type SNP --exclude-filtered -O {output[1]}"
