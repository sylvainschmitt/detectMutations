rule cp_reference:
    input:
        expand("{genome_dir}/{reference}.fa", genome_dir=config["genome_dir"], 
                allow_missing=True)
    output:
        "results/reference/{reference}.fa"
    log:
        "results/logs/cp_reference_{reference}.log"
    benchmark:
        "results/benchmarks/cp_reference_{reference}.benchmark.txt"
    shell:
        "cp {input} {output}"
