rule strelka2tsv:
    input:
        "results/mutations/B{branch}_T{tip}_on_{chromosome}.tip.vcf"
    output:
        "results/mutations/B{branch}_T{tip}_on_{chromosome}.tip.tsv"
    log:
        "results/logs/strelka2tsv_B{branch}_T{tip}_on_{chromosome}.log"
    benchmark:
        "results/benchmarks/strelka2tsv_B{branch}_T{tip}_on_{chromosome}.benchmark.txt"
    script:
        "../scripts/strelka2tsv.R"
