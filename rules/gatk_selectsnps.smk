rule gatk_selectsnps:
    input:
        "results/mutations/{tumor}_vs_{normal}.gatk.raw.vcf"
    output:
        "results/mutations/{tumor}_vs_{normal}_gatk.snps.vcf"
    log:
        "results/logs/gatk_selectsnps_{tumor}_vs_{normal}.log"
    benchmark:
        "results/benchmarks/gatk_selectsnps_{tumor}_vs_{normal}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\"  SelectVariants -V {input} -select-type SNP -O {output}"
        