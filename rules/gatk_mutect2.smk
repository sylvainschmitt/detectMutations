base=[config["base"]]
mutated=[config["mutated"]]

rule gatk_mutect2:
    input:
        "results/alignments/bwa/{base}.md.bam",
        "results/alignments/bwa/{mutated}.md.bam",
         multiext("results/raw_data/reference/reference", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict")
    output:
        "results/somatic/gatk/{mutated}_vs_{base}.vcf"
    log:
        "results/logs/gatk_mutect2_{mutated}_vs_{base}.log"
    benchmark:
        "results/benchmarks/gatk_mutect2_{mutated}_vs_{base}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk Mutect2 -R {input[2]} -I {input[1]} -tumor {wildcards.mutated}  -I {input[0]} -normal {wildcards.base}"
        " -dont-use-soft-clipped-bases true -O {output}"
