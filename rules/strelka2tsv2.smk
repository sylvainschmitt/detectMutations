rule strelka2tsv2:
    input:
        "results/mutations_inter/{tumor}_vs_{normal}_on_{chromosome}_{caller}.vcf"
    output:
        "results/mutations_tsv2/{tumor}_vs_{normal}_on_{chromosome}_{caller}.tsv"
    log:
        "results/logs/strelka2tsv_{tumor}_{normal}_{chromosome}_{caller}.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{tumor}_{normal}_{chromosome}_{caller}.benchmark.txt"
    script:
        "../scripts/strelka2tsv.R"
