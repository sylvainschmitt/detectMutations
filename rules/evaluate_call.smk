rule evaluate_call:
    input:
        "results/mutations/{lib}_AF{AF}_NR{NR}_REP{REP}_{caller}.vcf",
        expand("{mutdir}{reference}_REP{REP}_mutated_{lib}.tsv", mutdir=config["mutdir"], reference=config["reference"], allow_missing=True)
    output:
        "results/stats/{lib}_AF{AF}_NR{NR}_REP{REP}_{caller}.tsv"
    log:
        "results/logs/evaluate_call_{lib}_AF{AF}_NR{NR}_REP{REP}_{caller}.log"
    benchmark:
        "results/benchmarks/evaluate_call_{lib}_AF{AF}_NR{NR}_REP{REP}_{caller}.benchmark.txt"
    singularity:
        "https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif"
    script:
        "../scripts/evaluate_call.R"
