rule samtools_sort:
    input:
        "results/{library}/{library}.sam"
    output:
        temp("results/{library}/{library}.bam")
    log:
        "results/logs/samtools_sort_{library}.log"
    benchmark:
        "results/benchmarks/samtools_sort_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools sort --threads {threads} {input} > {output}"
