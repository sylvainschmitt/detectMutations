rule varscan2vcf:
    input:
        "results/{lib}_REP{REP}/varscan/{lib}_REP{REP}.snp",
        expand("{refdir}{reference}_REP{REP}_snps.vcf", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    output:
        temp("results/{lib}_REP{REP,\d+}/varscan/{lib}_REP{REP}.unfiltered.vcf")
    log:
        "results/logs/varscan2vcf_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/varscan2vcf_{lib}_REP{REP}.benchmark.txt"
    singularity:
        "https://github.com/sylvainschmitt/singularity-tidyverse-Biostrings/releases/download/0.0.1/sylvainschmitt-singularity-tidyverse-Biostrings.latest.sif"
    script:
        "../scripts/varscan2vcf.R"
