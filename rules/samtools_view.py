rule samtools_view:
    input:
        "results/alns/{library}_on_{reference}.sam",
        "results/reference/{reference}.fa"
    output:
        temp("results/alns/{library}_on_{reference}.raw.cram")
    log:
        "results/logs/samtools_view_{library}_on_{reference}.log"
    benchmark:
        "results/benchmarks/samtools_view_{library}_on_{reference}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools view -C -T {input[1]} -f 1 -F 12 {input[0]} > {output}"