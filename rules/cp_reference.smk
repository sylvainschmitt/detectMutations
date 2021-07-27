rule cp_reference:
    input:
        expand("{refdir}/{reference}.fa", refdir=config["refdir"], reference=config["reference"])
    output:
        expand("results/reference/{reference}.fa", reference=config["reference"])
    log:
        "results/logs/cp_reference.log"
    benchmark:
        "results/benchmarks/cp_reference.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "cp {input} {output}"
