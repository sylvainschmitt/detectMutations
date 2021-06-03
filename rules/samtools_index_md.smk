rule samtools_index_md:
    input:
        "results/{library}/{library}_{chromosome}.md.cram"
    output:
        "results/{library}/{library}_{chromosome}.md.cram.crai"
    log:
        "results/logs/samtools_index_md_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/samtools_index_md_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input}"
