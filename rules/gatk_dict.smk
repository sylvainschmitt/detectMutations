rule gatk_dict:
    input:
        "results/raw_data/reference/reference.fa"
    output:
        "results/raw_data/reference/reference.dict"
    log:
        "results/logs/gatk_dict.log"
    benchmark:
        "results/benchmarks/gatk_dict.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk CreateSequenceDictionary R={input} O={output}"
