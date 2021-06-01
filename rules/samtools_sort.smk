rule samtools_sort:
    input:
        "results/{library}/{library}_{chromosome}.sam"
    output:
        temp("results/{library}/aln/{library}_{chromosome}.bam")
    log:
        "results/logs/samtools_sort_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/samtools_sort_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools sort --threads {threads} {input} > {output}"
