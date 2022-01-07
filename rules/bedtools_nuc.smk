rule bedtools_nuc:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/reference/{reference}.1kb.bed", reference=config["reference"])
    output:
        expand("results/reference/{reference}.gc", reference=config["reference"])
    log:
        "results/logs/bedtools_nuc.log"
    benchmark:
        "results/benchmarks/bedtools_nuc.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools nuc -fi {input[0]} -bed  {input[1]} > {output}"
        