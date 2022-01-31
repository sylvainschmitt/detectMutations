rule strelka2tsv:
    input:
        "results/mutations/{tumor}_vs_{base}_on_{reference}.raw.vcf"
    output:
        "results/mutations/{tumor}_vs_{base}_on_{reference}.raw.tsv"
    log:
        "results/logs/strelka2tsv_{tumor}_vs_{base}_{reference}.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{tumor}_vs_{base}_{reference}.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    params:
        tumor="{tumor}",
        normal="{base}",
        reference="{reference}"
    script:
        "../scripts/strelka2tsv.R"
