rule filter_mutations:
    input:
        "results/mutations_{tissue}/{file}.sql"
    output:
        "results/mutations_{tissue}/{file}_filtered.tsv",
    log:
        "results/logs/filter_mutations_{tissue}_{file}.log"
    benchmark:
        "results/benchmarks/filter_mutations_{tissue}_{file}.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    threads: 1
    resources:
        mem_mb=1000
    params:
        lowDP=50,
        highDP=360,
        minAC=10,
        maxAF=0.8
    script:
        "../scripts/filter_mutations.R"
