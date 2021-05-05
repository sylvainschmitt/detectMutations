base=[config["base"]]
mutated=[config["mutated"]]

rule freebayes_somatic:
    input:
        "results/alignments/bwa/{base}.md.bam",
        "results/alignments/bwa/{mutated}.md.bam",
        multiext("results/raw_data/reference/reference", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict")
    output:
        "results/somatic/freebayes/{mutated}_vs_{base}.vcf"
    log:
        "results/logs/freebayes_somatic_{mutated}_vs_{base}.log"
    benchmark:
        "results/benchmarks/freebayes_somatic_{mutated}_vs_{base}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/freebayes/freebayes:latest"
    shell:
        "freebayes -f {input[2]} --pooled-continuous --pooled-discrete --genotype-qualities --report-genotype-likelihood-max"
        " --allele-balance-priors-off --min-alternate-fraction 0.03 --min-repeat-entropy 1 --min-alternate-count 2"
        " {input[0]} {input[1]} > {output}"
