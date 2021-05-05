sample=[config["base"], config["mutated"]]

rule fastqc:
    input:
         expand("results/raw_data/reads/{sample}_{strand}.fastq", strand=["R1","R2"], allow_missing=True)
    output:
        temp(expand("results/qc/reads/{sample}_{strand}_fastqc.{ext}", strand=["R1","R2"], ext=["html", "zip"], allow_missing=True))
    log:
        "results/logs/fastqc_{sample}.log"
    benchmark:
        "results/benchmarks/fastqc_{sample}.benchmark.txt"
    singularity: 
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    shell:
        "fastqc -t {threads} -q {input} ; "
        "mv results/raw_data/reads/{wildcards.sample}_R1_fastqc.html results/qc/reads/ ; "
        "mv results/raw_data/reads/{wildcards.sample}_R2_fastqc.html results/qc/reads/ ; "
        "mv results/raw_data/reads/{wildcards.sample}_R1_fastqc.zip results/qc/reads/ ; "
        "mv results/raw_data/reads/{wildcards.sample}_R2_fastqc.zip results/qc/reads/"
