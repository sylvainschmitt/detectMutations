rule bedtools_getfasta:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/mutations/{reference}_mutations.bed", reference=config["reference"])
    output:
        temp(expand("results/mutations/{reference}_mutations.fa", reference=config["reference"]))
    log:
        "results/logs/bedtools_getfasta.log"
    benchmark:
        "results/benchmarks/bedtools_getfasta.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools getfasta -fi {input[0]} -bed {input[1]} -fo {output}"