rule gatk_gathervcfs:
    input:
        expand("results/mutations/{tumor}_vs_{normal}_on_{interval}_gatk.raw.vcf", interval=intervals, allow_missing=True)
    output:
        "results/mutations/{tumor}_vs_{normal}.gatk.raw.vcf"
    log:
        "results/logs/gatk_gathervcfs_{tumor}_vs_{normal}.log"
    benchmark:
        "results/benchmarks/gatk_gathervcfs_{tumor}_vs_{normal}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        prefix=lambda wildcards, input: "-I " + " -I ".join(input),
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk GatherVcfs --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" {params.prefix} -O {output}"
        