rule manta:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{snps}", snps=config["snps"]),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
        expand("results/reference/{snps}.idx", snps=config["snps"])
    output:
        "results/{library}/manta/{library}.vcf",
        "results/{library}/manta/results/variants/candidateSmallIndels.vcf.gz"
    log:
        "results/logs/manta_{library}.log"
    benchmark:
        "results/benchmarks/manta_{library}.benchmark.txt"
    singularity: 
        "docker://quay.io/wtsicgp/strelka2-manta"
    shell:
        "rm -r results/{wildcards.library}/manta ; "
        "configManta.py "
        "--normalBam {input[2]} "
        "--tumorBam {input[1]} "
        "--referenceFasta {input[0]} "
        "--runDir results/{wildcards.library}/manta/ ; "
        "results/{wildcards.library}/manta/runWorkflow.py -m local -j {threads} ; "
        "zcat results/{wildcards.library}/manta/results/variants/somaticSV.vcf.gz > {output}"
        