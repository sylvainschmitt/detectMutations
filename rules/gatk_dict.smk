rule gatk_dict:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    output:
        expand("{refdir}{reference}_REP{REP,\d+}.dict", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    log:
        "results/logs/gatk_dict_REP{REP}.log"
    benchmark:
        "results/benchmarks/gatk_dict_REP{REP}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk CreateSequenceDictionary R={input} O={output}"
