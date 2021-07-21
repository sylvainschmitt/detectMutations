rule samtools_sort:
    input:
        "results/alns/{library}_on_{chromosome}.raw.cram"
    output:
        temp("results/alns/{library}_on_{chromosome}.sorted.cram")
    log:
        "results/logs/samtools_sort_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/samtools_sort_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools sort --threads {threads} -O cram {input} > {output}"
