rule freebayes:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        temp("results/{library}/freebayes/{library}.unfiltered.vcf")
    log:
        "results/logs/freebayes_{library}.log"
    benchmark:
        "results/benchmarks/freebayes_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/freebayes/freebayes:latest"
    shell:
        "freebayes -f {input[0]} --pooled-continuous --pooled-discrete --genotype-qualities --report-genotype-likelihood-max"
        " --allele-balance-priors-off --min-alternate-fraction 0.03 --min-repeat-entropy 1 --min-alternate-count 2"
        " {input[2]} {input[1]} > {output}"
