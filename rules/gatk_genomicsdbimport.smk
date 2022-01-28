rule gatk_genomicsdbimport:
    input:
        expand("results/hz/{library}.g.vcf", library=config["hz"])
    output:
        directory("results/hz/db")
    log:
        "results/logs/gatk_genomicsdbimport.log"
    benchmark:
        "results/benchmarks/gatk_genomicsdbimport.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        gvcfs=lambda wildcards, input: "-V " + " -V ".join(sorted(input)),
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk GenomicsDBImport --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" "
        "{params.gvcfs} --genomicsdb-workspace-path {output}"
