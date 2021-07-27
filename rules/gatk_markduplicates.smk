rule gatk_markduplicates:
    input:
        "results/{library}/{library}.sorted.cram",
         expand("results/reference/{reference}.fa", reference=config["reference"]),
        "results/{library}/{library}.sorted.cram.crai"
    output:
        temp("results/{library}/{library}.md.bam"),
        temp("results/{library}/{library}.md.bam.metrics")
    log:
        "results/logs/gatk_markduplicates_{library}.log"
    benchmark:
        "results/benchmarks/gatk_markduplicates_{library}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=100000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk MarkDuplicates --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -I {input[0]} -O {output[0]} -M {output[1]} -R {input[1]}"
