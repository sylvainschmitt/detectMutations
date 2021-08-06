rule strelka2tsv:
    input:
        "results/mutations/{tumor}_vs_{normal}_strelka2.vcf"
    output:
        temp("results/mutations/{tumor}_vs_{normal}_strelka2.tsv")
    log:
        "results/logs/strelka2tsv_{tumor}_{normal}_strelka2.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{tumor}_{normal}_strelka2.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.2/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/strelka2tsv.R"
        