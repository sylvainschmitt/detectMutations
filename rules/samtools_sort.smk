rule samtools_sort:
    input:
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.sam", type=["mutated", "base"], allow_missing=True)
    output:
        temp(expand("results/{lib}_REP{REP,\d+}/{lib}_REP{REP}_{type}.bam", type=["mutated", "base"], allow_missing=True))
    log:
        "results/logs/samtools_sort_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/samtools_sort_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools sort --threads {threads} {input[0]} > {output[0]} ; "
        "samtools sort --threads {threads} {input[1]} > {output[1]}"
