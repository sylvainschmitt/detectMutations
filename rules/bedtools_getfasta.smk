rule bedtools_getfasta:
    input:
        expand("{libdir}/haplome_v2.3.fa", libdir=config["libdir"]),
        "results/3P/mutations.bed"
    output:
        temp("results/3P/mutations.fa")
    log:
        "results/logs/bedtools_getfasta.log"
    benchmark:
        "results/benchmarks/bedtools_getfasta.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools getfasta -fi {input[0]} -bed {input[1]} -fo {output}"
        