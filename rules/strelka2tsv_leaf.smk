rule strelka2tsv_leaf:
    input:
        expand("results/mutations/{leaf}_vs_{base}.nontrunk.vcf", base=config["base"], allow_missing=True)
    output:
        expand("results/mutations/{leaf}_vs_{base}.nontrunk.tsv", base=config["base"], allow_missing=True)
    log:
        "results/logs/strelka2tsv_{leaf}.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{leaf}.benchmark.txt"
    params:
        tumor="{leaf}",
        normal=config["base"]
    script:
        "../scripts/strelka2tsv.R"
