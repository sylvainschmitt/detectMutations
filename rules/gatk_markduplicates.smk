rule gatk_markduplicates:
    input:
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.bam.bai", type=["mutated", "base"], allow_missing=True)
    output:
        expand("results/{lib}_REP{REP,\d+}/{lib}_REP{REP}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        temp(expand("results/{lib}_REP{REP,\d+}/{lib}_REP{REP}_{type}.bam.metrics", type=["mutated", "base"], allow_missing=True))
    log:
        "results/logs/gatk_markduplicates_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/gatk_markduplicates_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    resources:
        mem_mb=16000
    shell:
        "gatk MarkDuplicates --java-options \"-Xmx16G -Xms1G -Djava.io.tmpdir=tmp\" I={input[0]} O={output[0]} M={output[2]} ; "
        "gatk MarkDuplicates --java-options \"-Xmx16G -Xms1G -Djava.io.tmpdir=tmp\" I={input[1]} O={output[1]} M={output[3]}"
