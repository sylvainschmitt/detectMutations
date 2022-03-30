rule strelka2:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        "results/alns/{tumor}.md.cram",
        "results/alns/{base}.md.cram",
        "results/alns/{tumor}.md.cram.crai",
        "results/alns/{base}.md.cram.crai"
    output:
        "results/mutations_{tissue}/{tumor}_vs_{base}.raw.vcf"
    log:
        "results/logs/strelka2_{tissue}_{tumor}_vs_{base}.log"
    benchmark:
        "results/benchmarks/strelka2_{tissue}_{tumor}_vs_{base}.benchmark.txt"
    singularity: 
        "docker://quay.io/wtsicgp/strelka2-manta"
    threads: 20
    resources:
        mem_mb=100000
    shell:
        "configureStrelkaSomaticWorkflow.py "
        "--normalBam {input[2]} "
        "--tumorBam {input[1]} "
        "--referenceFasta {input[0]} "
        "--runDir tmp/strelka2_{wildcards.tumor}_vs_{wildcards.base} ; "
        "tmp/strelka2_{wildcards.tumor}_vs_{wildcards.base}/runWorkflow.py -m local -j {threads} ; "
        "zcat tmp/strelka2_{wildcards.tumor}_vs_{wildcards.base}/results/variants/somatic.snvs.vcf.gz > {output} ; "
        "rm -r tmp/strelka2_{wildcards.tumor}_vs_{wildcards.base}"
