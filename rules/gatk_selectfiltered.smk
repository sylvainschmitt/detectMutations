rule gatk_selectfiltered:
    input:
        "results/mutations/{tumor}_vs_{normal}_gatk.filtered.vcf"
    output:
        "results/mutations/{tumor}_vs_{normal}_gatk.vcf"
    log:
        "results/logs/gatk_selectfiltered_{tumor}_vs_{normal}.log"
    benchmark:
        "results/benchmarks/gatk_selectfiltered_{tumor}_vs_{normal}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\"  SelectVariants -V {input} --exclude-filtered -O {output}"
        