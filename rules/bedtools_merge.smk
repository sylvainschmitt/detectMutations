rule bedtools_merge:
    input:
        expand("results/mutations/{trunk}_vs_{base}.raw.vcf", trunk=config["trunk"], base=config["base"])
    output:
        expand("results/mutations/Tall_vs_{base}.raw.vcf", base=config["base"])
    log:
        "results/logs/bedtools_merge_Tall.log"
    benchmark:
        "results/benchmarks/bedtools_merge_Tall.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools merge -header -i {input} > {output}"
