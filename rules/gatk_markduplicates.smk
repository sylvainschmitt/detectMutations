rule gatk_markduplicates:
    input:
        "results/{library}/{library}_{chromosome}.sorted.cram",
         expand("results/reference/{reference}_{chromosome}.fa", reference=config["reference"], allow_missing = True),
        "results/{library}/{library}_{chromosome}.sorted.cram.crai"
    output:
        temp("results/{library}/{library}_{chromosome}.md.bam"),
        temp("results/{library}/{library}_{chromosome}.md.bam.metrics")
    log:
        "results/logs/gatk_markduplicates_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/gatk_markduplicates_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk MarkDuplicates --java-options "-Xmx10G -Xms100G -Djava.io.tmpdir=tmp" -I {input[0]} -O {output[0]} -M {output[1]} -R {input[1]}"
