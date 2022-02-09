rule strelka2sql_cambium:
    input:
        expand("results/mutations_cambium/{comps}.nonhz.tsv", comps=config["cambium_comp"])
    output:
        temp("results/mutations_cambium/cambium_nonhz_mutations.csv"),
        "results/mutations_cambium/cambium_nonhz_mutations.sql"
    log:
        "results/logs/strelka2sql_cambium.log"
    benchmark:
        "results/benchmarks/strelka2sql_cambium.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/strelka2sql.R"
