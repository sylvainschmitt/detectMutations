rule gatk_genotypegvcfs:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        directory("results/hz/db")
    output:
        "results/hz/raw_hz.vcf"
    log:
        "results/logs/gatk_genotypegvcfs.log"
    benchmark:
        "results/benchmarks/gatk_genotypegvcfs.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk GenotypeGVCFs --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" "
        "-R {input[0]} -V gendb://{input[1]} -O {output}"
        