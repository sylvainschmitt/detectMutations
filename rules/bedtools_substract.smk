caller2filter=["freebayes", "gatk"]

rule bedtools_substract:
    input:
        "results/{library}/{caller2filter}/{library}.unfiltered.vcf",
        expand("results/reference/{snps}", snps=config["snps"])
    output:
        "results/{library}/{caller2filter}/{library}.vcf"
    log:
        "results/logs/bedtools_substract_{library}_{caller2filter}.log"
    benchmark:
        "results/benchmarks/bedtools_substract_{library}_{caller2filter}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools subtract -header -a {input[0]} -b {input[1]} > {output}"
