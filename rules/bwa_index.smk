rule bwa_index:
    input:
        expand("results/reference/{reference}_{chromosome}.fa", reference=config["reference"], allow_missing = True)
    output:
        expand("results/reference/{reference}_{chromosome}.fa{ext}", 
                reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing = True)
    log:
        "results/logs/bwa_index_{chromosome}.log"
    benchmark:
        "results/benchmarks/bwa_index_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    shell:
        "bwa index {input}"
