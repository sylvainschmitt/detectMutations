rule mutect2:
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
        "results/mutations/{tumor}_vs_{normal}_mutect2.vcf"
    log:
        "results/logs/mutect2_{tumor}_vs_{normal}.log"
    benchmark:
        "results/benchmarks/mutect2_{tumor}_vs_{normal}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=100000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk Mutect2 --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\""
        "-R {input[0]} -I {input[1]} -tumor {wildcards.tumor}  -I {input[2]} -normal {wildcards.normal} "
        "-dont-use-soft-clipped-bases true -O {output}"
        