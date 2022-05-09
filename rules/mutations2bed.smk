rule mutations2bed:
    input:
        "results/mutations/mutations_filtered.tsv",
        expand("results/reference/{reference}.fa.fai", reference=config["reference"])
    output:
        temp(expand("results/mutations/{reference}_mutations.bed", reference=config["reference"])),
        temp(expand("results/mutations/{reference}_mutations.tsv", reference=config["reference"]))
    log:
        "results/logs/mutations2bed.log"
    benchmark:
        "results/benchmarks/mutations2bed.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    params:
        ref=config["reference"],
        N=100
    script:
        "../scripts/mutations2bed.R"
