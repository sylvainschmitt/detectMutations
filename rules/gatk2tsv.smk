rule gatk2tsv:
    input:
        expand("results/mutations/{vcfs}_gatk.vcf", vcfs=config["vcfs"])
    output:
        "results/mutations/gatk.raw.tsv"
    log:
        "results/logs/gatk2tsv.log"
    benchmark:
        "results/benchmarks/gatk2tsv.benchmark.txt"
    script:
        "../scripts/gatk2tsv.R"
