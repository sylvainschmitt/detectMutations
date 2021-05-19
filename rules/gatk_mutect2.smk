rule gatk_mutect2:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        "results/{library}/{library}.md.bam",
        expand("results/reference/{snps}", snps=config["snps"]),
        "results/{library}/{library}.md.bam.bai",
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
        expand("results/reference/{snps}.idx", snps=config["snps"])
    output:
        "results/{library}/mutect2/{library}.vcf"
    log:
        "results/logs/gatk_mutect2_{library}.log"
    benchmark:
        "results/benchmarks/gatk_mutect2_{library}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk Mutect2 -R {input[0]} -I {input[1]} --panel-of-normals {input[2]} -O {output}"

# "gatk Mutect2 -R {input[2]} -I {input[1]} -tumor {wildcards.mutated}  -I {input[0]} -normal {wildcards.base} -dont-use-soft-clipped-bases true -O {output}"
