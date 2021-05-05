sample=[config["base"], config["mutated"]]

rule qualimap:
    input:
        "results/alignments/bwa/{sample}.md.bam"
    output:
        "results/alignments/bwa/{sample}/qualimapReport.html"
    log:
        "results/logs/qualimap_{sample}.log"
    benchmark:
        "results/benchmarks/qualimap_{sample}.benchmark.txt"
    singularity: 
        "docker://pegi3s/qualimap"
    shell:
        "qualimap bamqc -bam {input} --paint-chromosome-limits -nt {threads} -skip-duplicated --skip-dup-mode 0 -outdir results/alignments/bwa/{wildcards.sample}  -outformat HTML"
