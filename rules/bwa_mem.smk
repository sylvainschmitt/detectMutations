rule bwa_mem:
    input:
        expand("results/reference/{reference}_{chromosome}.fa", reference=config["reference"], allow_missing = True),
        expand("results/{library}/{library}_{strand}.trimmed.paired.fastq.gz", strand=["1", "2"], allow_missing=True),
        expand("results/reference/{reference}_{chromosome}{ext}", 
                reference=config["reference"], ext=[".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict"],
                allow_missing=True)
    output:
        temp("results/{library}/{library}_{chromosome}.sam")
    params:
        rg=r"@RG\tID:{library}_{chromosome}\tSM:{library}_{chromosome}"
    log:
        "results/logs/bwa_mem_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/bwa_mem_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    threads: 8
    resources:
        mem_mb=32000
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input[0]} {input[1]} {input[2]} > {output}"
        