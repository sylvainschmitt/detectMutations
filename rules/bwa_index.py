rule bwa_index:
    input:
        "results/reference/{reference}.fa"
    output:
        expand("results/reference/{reference}.fa{ext}", allow_missing=True,
                ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    log:
        "results/logs/bwa_index_{reference}.log"
    benchmark:
        "results/benchmarks/bwa_index_{reference}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    shell:
        "bwa index {input}"