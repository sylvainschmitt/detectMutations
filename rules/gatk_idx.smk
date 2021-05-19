rule gatk_idx:
    input:
        expand("results/reference/{snps}", snps=config["snps"])
    output:
         expand("results/reference/{snps}.idx", snps=config["snps"])
    log:
        "results/logs/gatk_idx.log"
    benchmark:
        "results/benchmarks/gatk_idx.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk IndexFeatureFile -I {input}"
