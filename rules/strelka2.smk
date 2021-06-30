rule strelka2:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        temp("results/{library}/strelka2/{library}.unfiltered.vcf")
    log:
        "results/logs/strelka2_{library}.log"
    benchmark:
        "results/benchmarks/strelka2_{library}.benchmark.txt"
    singularity: 
        "docker://quay.io/wtsicgp/strelka2-manta"
    threads: 4
    resources:
        mem_mb=16000
    shell:
        "rm -r results/{wildcards.library}/strelka2 ; "
        "configureStrelkaSomaticWorkflow.py "
        "--normalBam {input[2]} "
        "--tumorBam {input[1]} "
        "--referenceFasta {input[0]} "
        "--runDir results/{wildcards.library}/strelka2/ ; "
        "results/{wildcards.library}/strelka2/runWorkflow.py -m local -j {threads} ; "
        "zcat results/{wildcards.library}/strelka2/results/variants/somatic.snvs.vcf.gz > {output}"
        
