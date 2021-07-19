rule gatk2tsv:
    input:
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_gatk.vcf"
    output:
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_gatk.tsv"
    log:
        "results/logs/gatk2tsv_{tumor}_{normal}_{chromosome}.log"
    benchmark:
        "results/benchmarks/gatk2tsv_{tumor}_{normal}_{chromosome}.benchmark.txt"
    script:
        "../scripts/gatk2tsv.R"
