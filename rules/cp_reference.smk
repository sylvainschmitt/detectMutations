rule cp_reference:
    input:
        expand("{refdir}{reference}.fa", refdir=config["refdir"], allow_missing=True)
    output:
        "results/reference/{reference}.fa"
    log:
        "results/logs/cp_{reference}.log"
    benchmark:
        "results/benchmarks/cp_{reference}.benchmark.txt"
    shell:
        "cp {input} {output}"
