rule gatk_haplotypecaller:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        "results/alns/{library}.md.cram",
        expand("{intervals}/{interval}", intervals=config["intervals"], allow_missing=True),
        expand("results/reference/{reference}.dict", reference=config["reference"])
    output:
        temp("results/hz/{library}/{interval}.g.vcf")
    log:
        "results/logs/gatk_haplotypecaller_{library}_{interval}.log"
    benchmark:
        "results/benchmarks/gatk_haplotypecaller_{library}_{interval}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk HaplotypeCaller --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" "
        "-R {input[0]} -I {input[1]} -L {input[2]} -O {output} -ERC GVCF"
