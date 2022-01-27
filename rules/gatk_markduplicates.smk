rule gatk_markduplicates:
    input:
        "results/alns/{library}_on_{reference}.sorted.cram",
         "results/reference/{reference}.fa",
        "results/alns/{library}_on_{reference}.sorted.cram.crai"
    output:
        temp("results/alns/{library}_on_{reference}.md.bam"),
        temp("results/alns/{library}_on_{reference}.md.bam.metrics")
    log:
        "results/logs/gatk_markduplicates_{library}_{reference}.log"
    benchmark:
        "results/benchmarks/gatk_markduplicates_{library}_{reference}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=100000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk MarkDuplicates --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -I {input[0]} -O {output[0]} -M {output[1]} -R {input[1]}"
