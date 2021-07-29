rule strelka2:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        "results/alns/B{branch}_T{tip}_L{repetition}.md.cram",
        "results/alns/C{repetition}.md.cram",
        "results/alns/B{branch}_T{tip}_L{repetition}.md.cram.crai",
        "results/alns/C{repetition}.md.cram.crai"
    output:
        temp("results/mutations/B{branch}_T{tip}_L{repetition}.raw.vcf",)
    log:
        "results/logs/strelka2_B{branch}_T{tip}_L{repetition}.log"
    benchmark:
        "results/benchmarks/strelka2_B{branch}_T{tip}_L{repetition}.benchmark.txt"
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
        "--runDir results/strelka2_B{wildcards.branch}_T{wildcards.tip}_L{wildcards.repetition} ; "
        "results/strelka2_B{wildcards.branch}_T{wildcards.tip}_L{wildcards.repetition}/runWorkflow.py -m local -j {threads} ; "
        "zcat results/strelka2_B{wildcards.branch}_T{wildcards.tip}_L{wildcards.repetition}/results/variants/somatic.snvs.vcf.gz > {output} ;"
        "rm -r results/strelka2_B{wildcards.branch}_T{wildcards.tip}_L{wildcards.repetition}"
