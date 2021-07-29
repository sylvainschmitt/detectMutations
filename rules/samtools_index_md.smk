rule samtools_index_md:
    input:
        "results/{library}/{library}.md.cram"
    output:
        "results/{library}/{library}.md.cram.crai"
    log:
        "results/logs/samtools_index_md_{library}.log"
    benchmark:
        "results/benchmarks/samtools_index_md_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input}"
