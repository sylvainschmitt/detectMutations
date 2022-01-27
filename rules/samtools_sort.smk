rule samtools_sort:
    input:
        "results/alns/{library}_on_{reference}.raw.cram"
    output:
        temp("results/alns/{library}_on_{reference}.sorted.cram")
    log:
        "results/logs/samtools_sort_{library}_{reference}.log"
    benchmark:
        "results/benchmarks/samtools_sort_{library}_{reference}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools sort --threads {threads} -O cram {input} > {output}"
