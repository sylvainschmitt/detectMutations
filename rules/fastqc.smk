rule fastqc:
    input:
        expand("{libdir}{library}_{type}_{strand}.raw.fastq", libdir=config["libdir"], type=["mutated", "base"], strand=["R1", "R2"], allow_missing=True)
    output:
        temp(expand("results/{library}/{library}_{type}_{strand}.raw_fastqc.{ext}", 
                    type=["mutated", "base"], strand=["R1", "R2"], ext=["html", "zip"], allow_missing=True))
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
        "fastqc -t {threads} -q {input} --outdir=results/{wildcards.library}/"
