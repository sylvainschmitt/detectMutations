rule samtools_view_md:
    input:
        "results/{library}/{library}.md.bam",
        expand("results/reference/{reference}.fa", reference=config["reference"], allow_missing = True)
    output:
        "results/{library}/{library}.md.cram"
    log:
        "results/logs/samtools_view_md_{library}.log"
    benchmark:
        "results/benchmarks/samtools_view_md_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools view -C -T {input[1]} {input[0]} > {output}"
