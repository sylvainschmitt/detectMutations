rule bcftools_shared:
    input:
        "results/hz/nonmissing_hz.vcf.gz"
    output:
        "results/hz/shared_hz.vcf.gz"
    log:
        "results/logs/bcftools_shared.log"
    benchmark:
        "results/benchmarks/bcftools_shared.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest"
    threads: 8
    resources:
        mem_mb=8000
    params:
        N=32
    shell:
        "bcftools view -i \'N_PASS(GT!=\"mis\" && GT!=\"RR\")>={params.N}\' --threads {threads} {input} | "
        "bgzip -c --threads {threads} > {output}"
