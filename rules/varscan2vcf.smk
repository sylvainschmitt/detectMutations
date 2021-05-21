rule varscan2vcf:
    input:
        "results/{library}/varscan/{library}.snp"
    output:
        "results/{library}/varscan/{library}.vcf"
    log:
        "results/logs/varscan2vcf_{library}.log"
    benchmark:
        "results/benchmarks/varscan2vcf_{library}.benchmark.txt"
    singularity:
        "https://github.com/sylvainschmitt/singularity-tidyverse-Biostrings/releases/download/0.0.1/sylvainschmitt-singularity-tidyverse-Biostrings.latest.sif"
    script:
        "../scripts/varscan2vcf.R"
