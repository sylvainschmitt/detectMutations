rule samtools_view_md:
    input:
        "results/alns/{library}_on_{chromosome}.md.bam",
        expand("results/reference/{reference}_{chromosome}.fa", reference=config["reference"], allow_missing = True)
    output:
        "results/alns/{library}_on_{chromosome}.md.cram"
    log:
        "results/logs/samtools_view_md_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/samtools_view_md_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools view -C -T {input[1]} {input[0]} > {output}"
