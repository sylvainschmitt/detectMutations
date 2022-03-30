rule samtools_index:
    input:
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.bam", type=["mutated", "base"], allow_missing=True)
    output:
        temp(expand("results/{lib}_REP{REP,\d+}/{lib}_REP{REP}_{type}.bam.bai", type=["mutated", "base"], allow_missing=True))
    log:
        "results/logs/samtools_index_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/samtools_index_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools index {input[0]} ; samtools index {input[1]}"
