rule samtools_faidx_napo:
    input:
        "results/napoleon/napoleon.fa"
    output:
        temp("results/napoleon/napoleon.fa.fai")
    log:
        "results/logs/samtools_faidx_napo.log"
    benchmark:
        "results/benchmarks/samtools_faidx_napo.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools faidx {input}"
