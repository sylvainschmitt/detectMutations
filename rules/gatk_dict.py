rule gatk_dict:
    input:
        "results/reference/{reference}.fa"
    output:
        "results/reference/{reference}.dict"
    log:
        "results/logs/gatk_dict_{reference}.log"
    benchmark:
        "results/benchmarks/gatk_dict_{reference}.benchmark.txt"
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
        
