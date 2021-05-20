rule gatk_mutect2:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{snps}", snps=config["snps"]),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
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
        "gatk Mutect2 -R {input[0]} -I {input[1]} -tumor {wildcards.library}_mutated  -I {input[2]} -normal {wildcards.library}_base "
        " --panel-of-normals {input[3]} -dont-use-soft-clipped-bases true -O {output}"
        
        ## working version ##
        # "gatk Mutect2 -R {input[0]} -I {input[1]} -I {input[2]} "
        # " --panel-of-normals {input[3]} -dont-use-soft-clipped-bases true -O {output}"
