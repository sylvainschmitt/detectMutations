rule fastqc:
    input:
         expand("results/{library}/{library}_{strand}.raw.fastq.gz", strand=["1", "2"], allow_missing=True)
    output:
        expand("results/{library}/{library}_{strand}.raw_fastqc.{ext}", strand=["1", "2"], ext=["html", "zip"], allow_missing=True)
    log:
        "results/logs/fastqc_{library}.log"
    benchmark:
        "results/benchmarks/fastqc_{library}.benchmark.txt"
    singularity: 
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads: 4
    resources:
        mem_mb=16000
    shell:
        "fastqc -t {threads} -q {input}"
