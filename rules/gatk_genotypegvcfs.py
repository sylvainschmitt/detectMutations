rule gatk_genotypegvcfs:
    input:
        "results/variants/db_{reference}",
        "results/reference/{reference}.fa"
    output:
        "results/{reference}.vcf"
    log:
        "results/logs/gatk_genotypegvcfs_{reference}.log"
    benchmark:
        "results/benchmarks/gatk_genotypegvcfs_{reference}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=2000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk GenotypeGVCFs --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" "
        "-R {input[1]} -V gendb://{input[0]} -O {output[0]}"
