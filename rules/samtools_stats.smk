rule samtools_stats:
    input:
        "results/{library}/{library}_{chromosome}.md.bam"
    output:
        temp("results/{library}/{library}_{chromosome}_md.bam.stats.out")
    log:
        "results/logs/samtools_stats_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/samtools_stats_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools stats --threads {threads} {input} > {output}"
