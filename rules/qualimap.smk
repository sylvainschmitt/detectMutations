rule qualimap:
    input:
        "results/{library}/{library}.md.bam"
    output:
        "results/{library}/qc/qualimap/qualimapReport.html"
    log:
        "results/logs/qualimap_{library}.log"
    benchmark:
        "results/benchmarks/qualimap_{library}.benchmark.txt"
    singularity: 
        "docker://pegi3s/qualimap"
    shell:
        "qualimap bamqc -bam {input} --paint-chromosome-limits -nt {threads} -skip-duplicated --skip-dup-mode " 
        "0 -outdir results/{wildcards.library}/qc/qualimap -outformat HTML"
