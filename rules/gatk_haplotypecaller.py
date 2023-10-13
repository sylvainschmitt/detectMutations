rule gatk_haplotypecaller:
    input:
        "results/alns/{library}_on_{reference}.sorted.cram",
        "results/reference/{reference}.fa",
        "results/alns/{library}_on_{reference}.sorted.cram.crai",
        "results/reference/{reference}.dict"
    output:
        temp("results/variants/{library}_on_{reference}.gvcf")
    log:
        "results/logs/gatk_haplotypecaller_{library}_on_{reference}.log"
    benchmark:
        "results/benchmarks/gatk_haplotypecaller_{library}_on_{reference}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=4000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk HaplotypeCaller --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\""
        " -R {input[1]} -I {input[0]} -O {output} -ERC GVCF"
