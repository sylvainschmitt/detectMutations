rule mv_data:
    input:
        "data/reference.fa",
        expand("data/{sample}_{strand}.fastq", sample=[config["base"], config["mutated"]], strand=["R1", "R2"])
    output:
        "results/raw_data/reference/reference.fa",
        expand("results/raw_data/reads/{sample}_{strand}.fastq", sample=[config["base"], config["mutated"]], strand=["R1", "R2"])
    log:
        "results/logs/mv_data.log"
    benchmark:
        "results/benchmarks/mv_data.benchmark.txt"
    shell:
        "cp {input[0]} {output[0]} ; "
        "cp {input[1]} {output[1]} ; "
        "cp {input[2]} {output[2]} ; "
        "cp {input[3]} {output[3]} ; "
        "cp {input[4]} {output[4]} "
        