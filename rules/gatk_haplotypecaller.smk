rule gatk_haplotypecaller:
    input:
        expand("results/reference/{reference}_{chromosome}.fa", 
                reference=config["reference"], allow_missing=True),
        "results/{tumor}/{tumor}_{chromosome}.md.cram",
        "results/{normal}/{normal}_{chromosome}.md.cram",
        "results/{tumor}/{tumor}_{chromosome}.md.cram.crai",
        "results/{normal}/{normal}_{chromosome}.md.cram.crai",
        expand("results/reference/{reference}_{chromosome}.fa{ext}", 
                reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        temp("results/mutations/{tumor}_vs_{normal}_on_{chromosome}_gatk.g.vcf"),
        temp("results/mutations/{tumor}_vs_{normal}_on_{chromosome}_gatk.g.vcf.idx")
    log:
        "results/logs/gatk_haplotypecaller_{tumor}_vs_{normal}_{chromosome}.log"
    benchmark:
        "results/benchmarks/gatk_haplotypecaller_{tumor}_vs_{normal}_{chromosome}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=100000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk HaplotypeCaller --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -R {input[0]} -I {input[1]} -O {output[0]} -ERC GVCF"
        