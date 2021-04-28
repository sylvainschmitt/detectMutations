rule index:
    input:
        "results/source/{file}"
    output:
        "results/source/{file}.fai"
    log:
        "results/logs/index_{file}.log"
    benchmark:
        "results/benchmarks/index_{file}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools faidx {input}"
