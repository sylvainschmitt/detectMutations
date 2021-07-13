rule strelka2tsv:
    input:
        "results/mutations_unique/{tumor}_vs_{normal}_on_{chromosome}_{caller}.vcf"
    output:
        "results/mutations_tsv/{tumor}_vs_{normal}_on_{chromosome}_{caller}.tsv"
    log:
        "results/logs/strelka2tsv_{tumor}_{normal}_{chromosome}_{caller}.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{tumor}_{normal}_{chromosome}_{caller}.benchmark.txt"
    script:
        "../scripts/strelka2tsv.R"
