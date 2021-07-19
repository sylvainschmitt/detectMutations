rule bedtools_subtract:
    input:
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_gatk.vcf",
        "results/mutations/{normal}_vs_{tumor}_on_{chromosome}_gatk.vcf"
    output:
        "results/mutations_filtered/{tumor}_vs_{normal}_on_{chromosome}_gatk.vcf"
    log:
        "results/logs/bedtools_substract_{tumor}_{normal}.log"
    benchmark:
        "results/benchmarks/bedtools_substract_{tumor}_{normal}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools subtract -header -a {input[0]} -b {input[1]} > {output}"
