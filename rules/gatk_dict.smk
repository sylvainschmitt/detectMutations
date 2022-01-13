rule gatk_dict:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"])
    output:
        expand("results/reference/{reference}.dict", reference=config["reference"])
    log:
        "results/logs/gatk_dict.log"
    benchmark:
        "results/benchmarks/gatk_dict.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk CreateSequenceDictionary --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" "
        "-R {input} -O {output}"
        