sample=[config["base"]]

rule freebayes_germline:
    input:
        "results/alignments/bwa/{sample}.md.bam",
         multiext("results/raw_data/reference/reference", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict")
    output:
        "results/germline/freebayes/{sample}.vcf"
    log:
        "results/logs/freebayes_germline_{sample}.log"
    benchmark:
        "results/benchmarks/freebayes_germline_{sample}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/freebayes/freebayes:latest"
    shell:
        "freebayes -f {input[1]} --min-alternate-fraction 0.1 --min-mapping-quality 1 {input[0]} > {output}"
        