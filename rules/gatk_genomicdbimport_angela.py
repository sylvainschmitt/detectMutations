rule gatk_genomicdbimport_angela:
    input:
        gvcfs=expand("results/variants/{library}_on_angela.gvcf",
                      library=config["angela_libs"], allow_missing=True),
        interval="data/angela.list"
    output:
        directory("results/variants/db_angela")
    log:
        "results/logs/gatk_genomicdbimport_angela.log"
    benchmark:
        "results/benchmarks/gatk_genomicdbimport_angela.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=2000
    params:
        gvcfs=lambda wildcards, input: "-V " + " -V ".join(sorted(input.gvcfs)),
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk GenomicsDBImport --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" "
        "{params.gvcfs} -L {input.interval} --genomicsdb-workspace-path {output}"
