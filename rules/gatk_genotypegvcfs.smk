rule gatk_genotypegvcfs:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        directory("results/hz/db"),
        expand("{intervals}/{interval}", intervals=config["intervals"], allow_missing=True)
    output:
        "results/hz/raw_hz/{interval}.vcf",
         temp("results/hz/raw_hz/{interval}.vcf.idx")
    log:
        "results/logs/gatk_genotypegvcfs_{interval}.log"
    benchmark:
        "results/benchmarks/gatk_genotypegvcfs_{interval}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk GenotypeGVCFs --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" "
        "-R {input[0]} -V gendb://{input[1]} -L {input[2]} -O {output[0]} --max-alternate-alleles 200 --max-genotype-count 2048"
        