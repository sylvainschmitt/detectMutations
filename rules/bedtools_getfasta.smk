rule bedtools_getfasta:
    input:
        "results/napoleon/napoleon.fa",
        "results/napoleon/napoleon_mutations.bed"
    output:
        temp("results/napoleon/napoleon_mutations.fa")
    log:
        "results/logs/bedtools_getfasta.log"
    benchmark:
        "results/benchmarks/bedtools_getfasta.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools getfasta -fi {input[0]} -bed {input[1]} -fo {output}"
