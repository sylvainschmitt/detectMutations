rule strelka2tsv_trunk:
    input:
        expand("results/mutations/{trunk}_vs_{base}.raw.vcf", base=config["base"], allow_missing=True)
    output:
        expand("results/mutations/{trunk}_vs_{base}.raw.tsv", base=config["base"], allow_missing=True)
    log:
        "results/logs/strelka2tsv_{trunk}.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{trunk}.benchmark.txt"
    params:
        tumor="{trunk}",
        normal=config["base"]
    script:
        "../scripts/strelka2tsv.R"
