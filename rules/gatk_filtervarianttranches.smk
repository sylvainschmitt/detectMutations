rule gatk_filtervarianttranches:
    input:
        "results/mutations/{tumor}_vs_{normal}_gatk.annotated.vcf"
    output:
        "results/mutations/{tumor}_vs_{normal}_gatk.vcf"
    log:
        "results/logs/gatk_filtervarianttranches_{tumor}_vs_{normal}.log"
    benchmark:
        "results/benchmarks/gatk_filtervarianttranches_{tumor}_vs_{normal}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk FilterVariantTranches --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -V {input} --info-key CNN_1D --snp-tranche 99.95 -O {output}"
        