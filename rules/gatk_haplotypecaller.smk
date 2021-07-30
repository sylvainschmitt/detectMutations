rule gatk_haplotypecaller:
    input:
        expand("results/reference/{reference}.fa", 
                reference=config["reference"], allow_missing=True),
        "results/{tumor}/{tumor}.md.cram",
         expand("{refdir}/intervals/{interval}", refdir=config["refdir"], allow_missing=True),
        "results/{normal}/{normal}.md.cram",
        "results/{tumor}/{tumor}.md.cram.crai",
        "results/{normal}/{normal}.md.cram.crai",
        expand("results/reference/{reference}.fa{ext}", 
                reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        "results/mutations/{tumor}_vs_{normal}_on_{interval}_gatk.raw.vcf"
    log:
        "results/logs/gatk_haplotypecaller_{tumor}_vs_{normal}_on_{interval}.log"
    benchmark:
        "results/benchmarks/gatk_haplotypecaller_{tumor}_vs_{normal}_on_{interval}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk HaplotypeCaller --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -R {input[0]} -I {input[1]} -L {input[2]} -O {output}"
        