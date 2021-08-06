rule mutect2tsv:
    input:
        "results/mutations/{tumor}_vs_{normal}_mutect2.vcf"
    output:
        temp("results/mutations/{tumor}_vs_{normal}_mutect2.tsv")
    log:
        "results/logs/strelka2tsv_{tumor}_{normal}_mutect2.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{tumor}_{normal}_mutect2.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.2/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/mutect2tsv.R"
        