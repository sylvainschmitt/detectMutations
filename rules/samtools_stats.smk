rule samtools_stats:
    input:
        "results/alns/{library}_on_{reference}.md.cram"
    output:
        "results/alns/{library}_on_{reference}.md.cram.stats"
    log:
        "results/logs/samtools_stats_{library}_{reference}.log"
    benchmark:
        "results/benchmarks/samtools_stats_{library}_{reference}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools stats --threads {threads} {input} > {output} ; "
