rule bwa_index:
    input:
        "results/raw_data/reference/reference.fa"
    output:
        multiext("results/raw_data/reference/reference.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "results/logs/bwa_index.log"
    benchmark:
        "results/benchmarks/bwa_index.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    shell:
        "bwa index {input}"
