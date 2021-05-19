rule gatk_markduplicates:
    input:
        multiext("results/{library}/{library}", ".bam", ".bam.bai")
    output:
        "results/{library}/{library}.md.bam",
        temp("results/{library}/{library}.bam.metrics")
    log:
        "results/logs/gatk_markduplicates_{library}.log"
    benchmark:
        "results/benchmarks/gatk_markduplicates_{library}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk MarkDuplicates I={input[0]} O={output[0]} M={output[1]}"
