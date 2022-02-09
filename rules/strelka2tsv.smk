rule strelka2tsv:
    input:
        "results/mutations_{tissue}/{tumor}_vs_{base}.nonhz.vcf"
    output:
        "results/mutations_{tissue}/{tumor}_vs_{base}.nonhz.tsv"
    log:
        "results/logs/strelka2tsv_{tissue}_{tumor}_vs_{base}.log"
    benchmark:
        "results/benchmarks/strelka2tsv_{tissue}_{tumor}_vs_{base}.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    params:
        tumor="{leaf}",
        normal="{base}"
    script:
        "../scripts/strelka2tsv.R"
