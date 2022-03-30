rule freebayes:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}.fa{ext}", 
               refdir=config["refdir"], reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        temp("results/{lib}_REP{REP,\d+}/freebayes/{lib}_REP{REP}.unfiltered.vcf")
    log:
        "results/logs/freebayes_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/freebayes_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/freebayes/freebayes:latest"
    shell:
        "freebayes -f {input[0]} --pooled-continuous --pooled-discrete --genotype-qualities --report-genotype-likelihood-max"
        " --allele-balance-priors-off --min-alternate-fraction 0.03 --min-repeat-entropy 1 --min-alternate-count 2"
        " {input[2]} {input[1]} > {output}"
