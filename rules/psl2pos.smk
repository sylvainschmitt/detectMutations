rule psl2pos:
    input:
        expand("{libdir}/3P_mutations.tsv", libdir=config["libdir"]),
        "results/3P/mutations.psl"
    output:
        "results/3P/mutations.tsv"
    log:
        "results/logs/psl2pos.log"
    benchmark:
        "results/benchmarks/psl2pos.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.2/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/psl2pos.R"
        