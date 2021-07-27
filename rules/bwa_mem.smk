rule bwa_mem:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{strand}.trimmed.paired.fastq.gz", strand=["1", "2"], allow_missing=True),
        expand("results/reference/{reference}{ext}", 
                reference=config["reference"], ext=[".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict"])
    output:
        temp("results/{library}/{library}.sam")
    params:
        rg=r"@RG\tID:{library}\tSM:{library}"
    log:
        "results/logs/bwa_mem_{library}.log"
    benchmark:
        "results/benchmarks/bwa_mem_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    threads: 20
    resources:
        mem_mb=200000
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input[0]} {input[1]} {input[2]} > {output}"
        