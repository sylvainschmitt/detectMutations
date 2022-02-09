rule mutations2bed:
    input:
        "results/mutations/mutations.tsv",
        expand("results/reference/{reference}.fa.fai", reference=config["reference"])
    output:
        temp(expand("results/mutations/{reference}_mutations.bed", reference=config["reference"]))
    log:
        "results/logs/mutations2bed.log"
    benchmark:
        "results/benchmarks/mutations2bed.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    params:
        reference=config["reference"],
        N=500
    script:
        "../scripts/mutations2bed.R"
