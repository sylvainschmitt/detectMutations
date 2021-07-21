rule bedtools_intersect:
    input:
        expand("results/mutations/B{branch}_T{tip}_L{rep}_on_{chromosome}.nonhz.vcf", rep=config["repetitions"], allow_missing=True)
    output:
        "results/mutations/B{branch}_T{tip}_on_{chromosome}.tip.vcf"
    log:
        "results/logs/bedtools_intersect_B{branch}_T{tip}_on_{chromosome}.log"
    benchmark:
        "results/benchmarks/bedtools_intersect_B{branch}_T{tip}_on_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools intersect -header -a {input[0]} -b {input} > {output}"
