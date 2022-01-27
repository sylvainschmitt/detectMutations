rule bwa_mem:
    input:
        "results/reference/{reference}.fa",
        expand("results/{library}/{library}_{strand}.trimmed.paired.fq.gz", strand=["1", "2"], allow_missing=True),
        expand("results/reference/{reference}.{ext}", 
               ext=["fa.amb", "fa.ann", "fa.bwt", "fa.pac", "fa.sa"], allow_missing=True)
    output:
        temp("results/alns/{library}_on_{reference}.sam")
    params:
        rg=r"@RG\tID:{library}\tSM:{library}"
    log:
        "results/logs/bwa_mem_{library}_{reference}.log"
    benchmark:
        "results/benchmarks/bwa_mem_{library}_{reference}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    threads: 40
    resources:
        mem_mb=400000
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input[0]} {input[1]} {input[2]} > {output}"
        