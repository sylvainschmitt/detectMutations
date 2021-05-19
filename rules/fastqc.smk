rule fastqc:
    input:
         expand("results/{library}/{library}_{strand}.raw.fastq", strand=["R1", "R2"], allow_missing=True)
    output:
        temp(expand("results/{library}/{library}_{strand}.raw_fastqc.{ext}", strand=["R1", "R2"], ext=["html", "zip"], allow_missing=True))
    log:
        "results/logs/fastqc_{library}.log"
    benchmark:
        "results/benchmarks/fastqc_{library}.benchmark.txt"
    singularity: 
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    shell:
        "fastqc -t {threads} -q {input}"
