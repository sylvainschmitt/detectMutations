rule strelka2tsv:
    input:
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_strelka2.vcf"
    output:
        temp("results/mutations/{tumor}_vs_{normal}_on_{chromosome}_strelka2.tsv")
    log:
        "results/logs/strelka2tsv_{tumor}_{normal}_{chromosome}_strelka2.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{tumor}_{normal}_{chromosome}_strelka2.benchmark.txt"
    script:
        "../scripts/strelka2tsv.R"
