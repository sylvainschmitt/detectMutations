rule bedtools_substract:
    input:
        "results/{library}/freebayes/{library}.unfiltered.vcf",
        expand("results/reference/{snps}", snps=config["snps"])
    output:
        "results/{library}/freebayes/{library}.vcf"
    log:
        "results/logs/bedtools_substract_{library}.log"
    benchmark:
        "results/benchmarks/bedtools_substract_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools subtract -header -a {input[0]} -b {input[1]} > {output}"
