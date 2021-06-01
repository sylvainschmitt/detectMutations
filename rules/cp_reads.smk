rule cp_reads:
    input:
        expand("{libdir}/{library}_{strand}.fastq.gz", libdir=config["libdir"], strand=["1", "2"], allow_missing=True)
    output:
        temp(expand("results/{library}/{library}_{strand}.raw.fastq", strand=["1", "2"], allow_missing=True))
    log:
        "results/logs/cp_{library}.log"
    benchmark:
        "results/benchmarks/cp_{library}.benchmark.txt"
    shell:
        "zcat {input[0]} > {output[0]} ; zcat {input[1]} > {output[1]}"
        