rule samtools_index:
    input:
        "results/alns/{library}_on_{reference}.sorted.cram"
    output:
        temp("results/alns/{library}_on_{reference}.sorted.cram.crai")
    log:
        "results/logs/samtools_index_{library}_{reference}.log"
    benchmark:
        "results/benchmarks/samtools_index_{library}_{reference}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input}"
