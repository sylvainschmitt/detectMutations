rule bedtools_substract:
    input:
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_{caller}.vcf",
        "results/mutations/{normal}_vs_{tumor}_on_{chromosome}_{caller}.vcf"
    output:
        "results/mutations_unique/{tumor}_vs_{normal}_on_{chromosome}_{caller}.vcf"
    log:
        "results/logs/bedtools_substract_{tumor}_{normal}_{chromosome}_{caller}.log"
    benchmark:
        "results/benchmarks/bedtools_substract_{tumor}_{normal}_{chromosome}_{caller}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools subtract -header -a {input[0]} -b {input[1]} > {output}"
