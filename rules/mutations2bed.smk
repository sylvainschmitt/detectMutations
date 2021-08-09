rule mutations2bed:
    input:
        expand("{libdir}/3P_mutations.tsv", libdir=config["libdir"])
    output:
        temp("results/3P/mutations.bed")
    log:
        "results/logs/mutations2bed.log"
    benchmark:
        "results/benchmarks/mutations2bed.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.2/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/mutations2bed.R"
        