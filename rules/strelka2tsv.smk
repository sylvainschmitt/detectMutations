rule strelka2tsv:
    input:
        "results/mutations/{tumor}_vs_{base}.raw.vcf"
    output:
        temp("results/mutations/{tumor}_vs_{base}.raw.tsv")
    log:
        "results/logs/strelka2tsv_{tumor}_vs_{base}.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{tumor}_vs_{base}.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    params:
        tumor="{tumor}",
        normal="{base}"
    script:
        "../scripts/strelka2tsv.R"
