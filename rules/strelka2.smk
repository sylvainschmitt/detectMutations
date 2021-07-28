rule strelka2:
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
        "results/mutations/{tumor}_vs_{normal}_strelka2.vcf"
    log:
        "results/logs/strelka2_{tumor}_vs_{normal}.log"
    benchmark:
        "results/benchmarks/strelka2_{tumor}_vs_{normal}.benchmark.txt"
    singularity: 
        "docker://quay.io/wtsicgp/strelka2-manta"
    threads: 20
    resources:
        mem_mb=200000
    shell:
        "configureStrelkaSomaticWorkflow.py "
        "--normalBam {input[2]} "
        "--tumorBam {input[1]} "
        "--referenceFasta {input[0]} "
        "--runDir results/strelka2/{wildcards.tumor}_vs_{wildcards.normal} ; "
        "results/strelka2/{wildcards.tumor}_vs_{wildcards.normal}/runWorkflow.py -m local -j {threads} ; "
        "zcat results/strelka2/{wildcards.tumor}_vs_{wildcards.normal}/results/variants/somatic.snvs.vcf.gz > {output} ;"
        "rm -r results/strelka2/{wildcards.tumor}_vs_{wildcards.normal}"
        
