rule cp_reads:
    input:
        expand("{libdir}{library}_{type}_{strand}.fastq", libdir=config["libdir"], type=["mutated", "base"], strand=["R1", "R2"], allow_missing=True)
    output:
        temp(expand("results/{library}/{library}_{type}_{strand}.raw.fastq", type=["mutated", "base"], strand=["R1", "R2"], allow_missing=True))
    log:
        "results/logs/cp_{library}.log"
    benchmark:
        "results/benchmarks/cp_{library}.benchmark.txt"
    shell:
        "cp {input[0]} {output[0]} ; cp {input[1]} {output[1]} ; cp {input[2]} {output[2]} ; cp {input[3]} {output[3]}"
        