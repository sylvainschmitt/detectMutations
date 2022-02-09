rule bcftools_biallelic:
    input:
        "results/hz/raw_hz.vcf"
    output:
        temp("results/hz/biallelic_hz.vcf.gz")
    log:
        "results/logs/bcftools_biallelic.log"
    benchmark:
        "results/benchmarks/bcftools_biallelic.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest"
    threads: 8
    resources:
        mem_mb=8000
    shell:
        "bcftools view --max-alleles 2 --threads {threads} {input} | "
        "bgzip -c --threads {threads} > {output}"
