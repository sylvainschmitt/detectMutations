rule strelka2:
    input:
        expand("results/reference/{reference}_{chromosome}.fa", 
                reference=config["reference"], allow_missing=True),
        "results/{tumor}/{tumor}_{chromosome}.md.cram",
        "results/{normal}/{normal}_{chromosome}.md.cram",
        "results/{tumor}/{tumor}_{chromosome}.md.cram.crai",
        "results/{normal}/{normal}_{chromosome}.md.cram.crai",
        expand("results/reference/{reference}_{chromosome}.fa{ext}", 
                reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_strelka2.vcf"
    log:
        "results/logs/strelka2_{tumor}_vs_{normal}_{chromosome}.log"
    benchmark:
        "results/benchmarks/strelka2_{tumor}_vs_{normal}_{chromosome}.benchmark.txt"
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
        "--runDir results/strelka2/{wildcards.tumor}_vs_{wildcards.normal}_on_{wildcards.chromosome} ; "
        "results/strelka2/{wildcards.tumor}_vs_{wildcards.normal}_on_{wildcards.chromosome}/runWorkflow.py -m local -j {threads} ; "
        "zcat results/strelka2/{wildcards.tumor}_vs_{wildcards.normal}_on_{wildcards.chromosome}/results/variants/somatic.snvs.vcf.gz > {output}"
        
