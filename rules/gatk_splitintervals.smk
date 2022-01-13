rule gatk_splitintervals:
    input:
        expand("results/reference/{reference}.{ext}", reference=config["reference"], ext=["fa", "fa.fai", "dict"])
    output:
        directory("results/hz/intervals")
    log:
        "results/logs/gatk_splitintervals.log"
    benchmark:
        "results/benchmarks/gatk_splitintervals.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        n=config["intervals"],
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk SplitIntervals --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -R {input[0]} --scater-count {params.n} -O {output}"
