rule samtools_faidx_haplome:
    input:
        expand("{libdir}/haplome_v2.3.fa", libdir=config["libdir"])
    output:
        expand("{libdir}/haplome_v2.3.fa.fai", libdir=config["libdir"])
    log:
        "results/logs/samtools_faidx_haplome.log"
    benchmark:
        "results/benchmarks/samtools_faidx_haplome.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools faidx {input}"
