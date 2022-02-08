rule plink_nonmissing:
    input:
        "results/hz/intermediate_hz.vcf.gz"
    output:
        "results/hz/nonmissing_hz.vcf.gz",
        temp("results/hz/nonmissing_hz.log"),
        temp("results/hz/nonmissing_hz.nosex")
    log:
        "results/logs/plink_nonmissing.log"
    benchmark:
        "results/benchmarks/plink_nonmissing.benchmark.txt"
    singularity: 
        "docker://gelog/plink"
    threads: 1
    resources:
        mem_mb=1000
    params:
        prefix = "results/hz/nonmissing_hz"
    shell:
        "plink --vcf {input} --allow-extra-chr --mind 1 --geno 0 --recode vcf-iid --out {params.prefix} ; "
        "gzip {params.prefix}.vcf"
        