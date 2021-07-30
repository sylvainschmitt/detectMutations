rule gatk_cnnscorevariants:
    input:
        expand("results/reference/{reference}.fa", 
                reference=config["reference"], allow_missing=True),
        "results/mutations/{tumor}_vs_{normal}.gatk.raw.vcf",
        expand("results/reference/{reference}.fa{ext}", 
                reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        temp("results/mutations/{tumor}_vs_{normal}_gatk.annotated.vcf")
    log:
        "results/logs/gatk_cnnscorevariants_{tumor}_vs_{normal}.log"
    benchmark:
        "results/benchmarks/gatk_cnnscorevariants_{tumor}_vs_{normal}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk CNNScoreVariants --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -R {input[0]} -V {input[1]} -O {output}"
        