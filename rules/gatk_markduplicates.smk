rule gatk_markduplicates:
    input:
        "results/{library}/aln/{library}_{chromosome}.bam",
        "results/{library}/aln/{library}_{chromosome}.bam.bai"
    output:
        "results/{library}/{library}_{chromosome}_md.bam",
        temp("results/{library}/{library}_{chromosome}.bam.metrics")
    log:
        "results/logs/gatk_markduplicates_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/gatk_markduplicates_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk MarkDuplicates I={input[0]} O={output[0]} M={output[1]}"
