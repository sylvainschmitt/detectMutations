rule strelka2:
    input:
        "results/reference/{reference}.fa",
        "results/alns/{tumor}_on_{reference}.md.cram",
        "results/alns/{base}_on_{reference}.md.cram",
        "results/alns/{tumor}_on_{reference}.md.cram.crai",
        "results/alns/{base}_on_{reference}.md.cram.crai"
    output:
        "results/mutations/{tumor}_vs_{base}_on_{reference}.raw.vcf"
    log:
        "results/logs/strelka2_{tumor}_vs_{base}_{reference}.log"
    benchmark:
        "results/benchmarks/strelka2_{tumor}_vs_{base}_{reference}.benchmark.txt"
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
        "--runDir tmp/strelka2_{wildcards.tumor}_vs_{wildcards.base} ; "
        "tmp/strelka2_{wildcards.tumor}_vs_{wildcards.base}/runWorkflow.py -m local -j {threads} ; "
        "zcat tmp/strelka2_{wildcards.tumor}_vs_{wildcards.base}/results/variants/somatic.snvs.vcf.gz > {output} ; "
        "rm -r tmp/strelka2_{wildcards.tumor}_vs_{wildcards.base}"
