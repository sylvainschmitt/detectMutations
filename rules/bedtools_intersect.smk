rule bedtools_intersect:
    input:
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_{caller}.vcf",
        "results/mutations/{normal}_vs_{tumor}_on_{chromosome}_{caller}.vcf"
    output:
        "results/mutations_inter/{tumor}_vs_{normal}_on_{chromosome}_{caller}.vcf"
    log:
        "results/logs/bedtools_intersect_{tumor}_{normal}_{chromosome}_{caller}.log"
    benchmark:
        "results/benchmarks/bedtools_intersect_{tumor}_{normal}_{chromosome}_{caller}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools intersect -header -a {input[0]} -b {input[1]} > {output}"
