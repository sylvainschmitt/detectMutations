rule gatk_snps:
    input:
        "results/hz/biallelic_hz.vcf.gz"
    output:
        temp("results/hz/biallelic_hz.vcf.gz.tbi"),
        temp("results/hz/intermediate_hz.vcf.gz"),
        temp("results/hz/filtered_hz.vcf.gz")
    log:
        "results/logs/gatk_snps.log"
    benchmark:
        "results/benchmarks/gatk_snps.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=1000
    shell:
        "gatk IndexFeatureFile -I {input} ;"
        "gatk VariantFiltration -V {input} -filter-expression \"QUAL < 30.0 || QD < 2.0 || FS > 60.0 || SOR > 3.0\" --filter-name \"FAIL\" "
        "-O {output[1]} ;"
        "gatk SelectVariants -V {output[1]} -select-type SNP --exclude-filtered -O {output[2]}"
        