rule gatk_haplotypecaller:
    input:
        expand("results/reference/{reference}.fa", 
                reference=config["reference"], allow_missing=True),
        "results/{tumor}/{tumor}.md.cram",
        "results/{normal}/{normal}.md.cram",
        "results/{tumor}/{tumor}.md.cram.crai",
        "results/{normal}/{normal}.md.cram.crai",
        expand("results/reference/{reference}.fa{ext}", 
                reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        "results/mutations/{tumor}_vs_{normal}_gatk.raw.vcf"
    log:
        "results/logs/gatk_haplotypecaller_{tumor}_vs_{normal}.log"
    benchmark:
        "results/benchmarks/gatk_haplotypecaller_{tumor}_vs_{normal}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 20
    resources:
        mem_mb=100000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk HaplotypeCaller --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -R {input[0]} -I {input[1]} -O {output}"
        