rule blat:
    input:
        "results/reference/{references}.fa",
        expand("results/mutations/{reference}_mutations.fa", reference=config["reference"])
    output:
        temp(expand("results/mutations/{reference}_mutations_on_{references}.psl", reference=config["reference"], allow_missing=True))
    log:
        "results/logs/blat_{references}.log"
    benchmark:
        "results/benchmarks/blat_{references}.benchmark.txt"
    singularity: 
        "docker://quay.io/biocontainers/ucsc-blat:377--ha8a8165_4"
    shell:
        "blat {input[0]} {input[1]} {output}"