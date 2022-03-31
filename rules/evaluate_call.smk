rule evaluate_call:
    input:
        "results/mutations/{lib}_AF{AF}_NR{NR}_REP{REP}_{caller}.vcf",
        expand("{mutdir}{reference}_REP{REP}_mutated_{lib}.tsv", mutdir=config["mutdir"], reference=config["reference"], allow_missing=True)
    output:
        temp("results/stats/{lib}_AF{AF}_NR{NR}_REP{REP}_{caller}.tsv")
    log:
        "results/logs/evaluate_call_{lib}_AF{AF}_NR{NR}_REP{REP}_{caller}.log"
    benchmark:
        "results/benchmarks/evaluate_call_{lib}_AF{AF}_NR{NR}_REP{REP}_{caller}.benchmark.txt"
    singularity:
        "https://github.com/sylvainschmitt/singularity-tidyverse-Biostrings/releases/download/0.0.1/sylvainschmitt-singularity-tidyverse-Biostrings.latest.sif"
    script:
        "../scripts/evaluate_call.R"
