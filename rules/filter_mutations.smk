rule filter_mutations:
    input:
        "results/mutations_unique/{tumor}_vs_{normal}_on_{chromosome}_{caller}.vcf"
    output:
        "results/mutations_filtered/{tumor}_vs_{normal}_on_{chromosome}_{caller}.tsv"
    log:
        "results/logs/filter_mutations_{tumor}_{normal}_{chromosome}_{caller}.log"
    benchmark:
        "results/benchmarks/filter_mutations_{tumor}_{normal}_{chromosome}_{caller}.benchmark.txt"
    script:
        "../scripts/filter_mutations.R"
