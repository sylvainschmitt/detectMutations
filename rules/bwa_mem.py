rule bwa_mem:
    input:
        "results/reference/{reference}.fa",
        expand("results/reads/{library}/{library}_{strand}.trimmed.paired.fq.gz", 
                strand=["R1", "R2"], allow_missing=True),
        expand("results/reference/{reference}.{ext}", allow_missing=True,
                ext=["fa.amb", "fa.ann", "fa.bwt", "fa.pac", "fa.sa"])
    output:
        temp("results/alns/{library}_on_{reference}.sam")
    params:
        rg=r"@RG\tID:{library}\tSM:{library}"
    log:
        "results/logs/bwa_mem_{library}_on_{reference}.log"
    benchmark:
        "results/benchmarks/bwa_mem_{library}_on_{reference}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    threads: 4
    resources:
        mem_mb=2000
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input[0]} {input[1]} {input[2]} > {output}"
        
