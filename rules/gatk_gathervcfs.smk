rule gatk_gathervcfs:
    input:
        expand("results/hz/{library}/{interval}.g.vcf", interval=intervals, allow_missing=True)
    output:
        "results/hz/{library}.g.vcf"
    log:
        "results/logs/gatk_gathervcfs_{library}.log"
    benchmark:
        "results/benchmarks/gatk_gathervcfs_{library}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        gvcfs=lambda wildcards, input: "-I " + " -I ".join(sorted(input)),
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk GatherVcfs --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" {params.gvcfs} -O {output}"
