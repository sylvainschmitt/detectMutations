rule gatk_dict:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"])
    output:
        expand("results/reference/{reference}.dict", reference=config["reference"])
    log:
        "results/logs/gatk_dict.log"
    benchmark:
        "results/benchmarks/gatk_dict.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk CreateSequenceDictionary R={input} O={output}"
