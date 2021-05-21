rule muse:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{snps}", snps=config["snps"]),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
        expand("results/reference/{snps}.idx", snps=config["snps"])
    output:
        "results/{library}/muse/{library}.vcf"
    log:
        "results/logs/muse_{library}.log"
    benchmark:
        "results/benchmarks/muse_{library}.benchmark.txt"
    singularity: 
        "docker://opengenomics/muse"
    shell:
        "muse.py -f {input[0]} --tumor-bam {input[1]} --normal-bam {input[2]} -O {output} -n {threads} -w results/{wildcards.library}/muse/"
        