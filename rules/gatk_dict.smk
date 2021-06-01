rule gatk_dict:
    input:
        expand("results/reference/{reference}_{chromosome}.fa", reference=config["reference"], allow_missing = True)
    output:
        expand("results/reference/{reference}_{chromosome}.dict", reference=config["reference"], allow_missing = True)
    log:
        "results/logs/gatk_dict_{chromosome}.log"
    benchmark:
        "results/benchmarks/gatk_dict_{chromosome}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk CreateSequenceDictionary R={input} O={output}"
