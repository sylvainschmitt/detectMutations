rule gatk2tsv:
    input:
        expand("results/mutations/{vcfs}_on_{chromosome}_gatk.vcf", vcfs=config["vcfs"], allow_missing = True)
    output:
        "results/mutations/{chromosome}_gatk.tsv"
    log:
        "results/logs/gatk2tsv_{chromosome}.log"
    benchmark:
        "results/benchmarks/gatk2tsv_{chromosome}.benchmark.txt"
    script:
        "../scripts/gatk2tsv.R"
