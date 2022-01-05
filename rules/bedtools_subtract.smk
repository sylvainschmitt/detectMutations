rule bedtools_subtract:
    input:
        expand("results/mutations/{leaf}_vs_{base}.raw.vcf", base=config["base"], allow_missing=True),
        expand("results/mutations/{trunk}_vs_{base}.raw.vcf", trunk=config["trunk"], base=config["base"])
    output:
        expand("results/mutations/{leaf}_vs_{base}.nontrunk.vcf", base=config["base"], allow_missing=True)
    log:
        "results/logs/bedtools_subtract_{leaf}.log"
    benchmark:
        "results/benchmarks/bedtools_subtract_{leaf}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools subtract -header -a {input[0]} -b {input[1]} |"
        "bedtools subtract -header -a stdin -b {input[2]} |"
        "bedtools subtract -header -a stdin -b {input[3]} >"
        "{output}"
        