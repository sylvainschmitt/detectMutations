rule mutect2tsv:
    input:
        "results/mutations/{tumor}_vs_{normal}_mutect2.vcf"
    output:
        temp("results/mutations/{tumor}_vs_{normal}_mutect2.tsv")
    log:
        "results/logs/strelka2tsv_{tumor}_{normal}_mutect2.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{tumor}_{normal}_mutect2.benchmark.txt"
    script:
        "../scripts/mutect2tsv.R"
        