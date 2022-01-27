rule samtools_view_md:
    input:
        "results/alns/{library}_on_{reference}.md.bam",
        "results/reference/{reference}.fa"
    output:
        "results/alns/{library}_on_{reference}.md.cram"
    log:
        "results/logs/samtools_view_md_{library}_{reference}.log"
    benchmark:
        "results/benchmarks/samtools_view_md_{library}_{reference}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools view -C -T {input[1]} {input[0]} > {output}"
