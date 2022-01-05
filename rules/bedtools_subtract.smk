rule bedtools_subtract:
    input:
        expand("results/mutations/{leaf}_vs_{base}.raw.vcf", base=config["base"], allow_missing=True),
        expand("results/mutations/Tall_vs_{base}.raw.vcf", base=config["base"])
    output:
        expand("results/mutations/{leaf}_vs_{base}.nontrunk.vcf", base=config["base"], allow_missing=True)
    log:
        "results/logs/bedtools_subtract_{leaf}.log"
    benchmark:
        "results/benchmarks/bedtools_subtract_{leaf}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools subtract -header -a {input[0]} -b {input[1]} > {output}"