rule qualimap:
    input:
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True)
    output:
        expand("results/{library}/qualimap_{type}/qualimapReport.html", type=["mutated", "base"], allow_missing=True)
    log:
        "results/logs/qualimap_{library}.log"
    benchmark:
        "results/benchmarks/qualimap_{library}.benchmark.txt"
    singularity: 
        "docker://pegi3s/qualimap"
    shell:
        "qualimap bamqc -bam {input[0]} --paint-chromosome-limits -nt {threads} -skip-duplicated --skip-dup-mode " 
        "0 -outdir results/{wildcards.library}/qualimap_mutated -outformat HTML ; "
        "qualimap bamqc -bam {input[1]} --paint-chromosome-limits -nt {threads} -skip-duplicated --skip-dup-mode " 
        "0 -outdir results/{wildcards.library}/qualimap_base -outformat HTML"
