rule strelka2:
    input:
        expand("results/reference/{reference}_{chromosome}.fa", 
                reference=config["reference"], allow_missing=True),
        "results/alns/B{branch}_T{tip}_L{repetition}_on_{chromosome}.md.cram",
        "results/alns/C{repetition}_on_{chromosome}.md.cram",
        "results/alns/B{branch}_T{tip}_L{repetition}_on_{chromosome}.md.cram.crai",
        "results/alns/C{repetition}_on_{chromosome}.md.cram.crai"
    output:
        temp("results/mutations/B{branch}_T{tip}_L{repetition}_on_{chromosome}.raw.vcf",)
    log:
        "results/logs/strelka2_B{branch}_T{tip}_L{repetition}_on_{chromosome}.log"
    benchmark:
        "results/benchmarks/strelka2_B{branch}_T{tip}_L{repetition}_on_{chromosome}.benchmark.txt"
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
        "--runDir results/strelka2_B{wildcards.branch}_T{wildcards.tip}_L{wildcards.repetition}_on_{wildcards.chromosome} ; "
        "results/strelka2_B{wildcards.branch}_T{wildcards.tip}_L{wildcards.repetition}_on_{wildcards.chromosome}/runWorkflow.py -m local -j {threads} ; "
        "zcat results/strelka2_B{wildcards.branch}_T{wildcards.tip}_L{wildcards.repetition}_on_{wildcards.chromosome}/results/variants/somatic.snvs.vcf.gz > {output} ;"
        "rm -r results/strelka2_B{wildcards.branch}_T{wildcards.tip}_L{wildcards.repetition}_on_{wildcards.chromosome}"
        
