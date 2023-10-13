rule gatk_genomicdbimport_sixto:
    input:
        gvcfs=expand("results/variants/{library}_on_sixto.gvcf",
                      library=config["sixto_libs"], allow_missing=True),
        interval="data/sixto.list"
    output:
        directory("results/variants/db_sixto")
    log:
        "results/logs/gatk_genomicdbimport_sixto.log"
    benchmark:
        "results/benchmarks/gatk_genomicdbimport_sixto.benchmark.txt"
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
