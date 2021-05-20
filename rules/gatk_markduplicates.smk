rule gatk_markduplicates:
    input:
        expand("results/{library}/{library}_{type}.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.bam.bai", type=["mutated", "base"], allow_missing=True)
    output:
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        temp(expand("results/{library}/{library}_{type}.bam.metrics", type=["mutated", "base"], allow_missing=True))
    log:
        "results/logs/gatk_markduplicates_{library}.log"
    benchmark:
        "results/benchmarks/gatk_markduplicates_{library}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk MarkDuplicates I={input[0]} O={output[0]} M={output[2]} ; "
        "gatk MarkDuplicates I={input[1]} O={output[1]} M={output[3]}"
