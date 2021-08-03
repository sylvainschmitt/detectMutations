rule gatk_mutect2:
    input:
        expand("results/reference/{reference}.fa", 
                reference=config["reference"], allow_missing=True),
        "results/{tumor}/{tumor}.md.cram",
        "results/{normal}/{normal}.md.cram",
        "results/reference/intervals/{interval}",
        "results/{tumor}/{tumor}.md.cram.crai",
        "results/{normal}/{normal}.md.cram.crai",
        expand("results/reference/{reference}.fa{ext}", 
                reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        temp("results/mutations/{tumor}_vs_{normal}_on_{interval}_mutect2.interval.vcf")
    log:
        "results/logs/gatk_mutect2_{tumor}_vs_{normal}_{interval}.log"
    benchmark:
        "results/benchmarks/gatk_mutect2_{tumor}_vs_{normal}_{interval}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk Mutect2 --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" "
        "-R {input[0]} -L {input[3]} -I {input[1]} -tumor {wildcards.tumor}  -I {input[2]} -normal {wildcards.normal} "
        "-dont-use-soft-clipped-bases true -O {output}"
        