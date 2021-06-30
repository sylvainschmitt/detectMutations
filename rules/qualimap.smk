rule qualimap:
    input:
        "results/{library}/{library}_{chromosome}.md.cram"
    output:
        "results/{library}/qualimap_{library}_{chromosome}/qualimapReport.html"
    log:
        "results/logs/qualimap_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/qualimap_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "docker://pegi3s/qualimap"
    threads: 4
    resources:
        mem_mb=16000
    shell:
        "qualimap bamqc -bam {input} --paint-chromosome-limits -nt {threads} -skip-duplicated --skip-dup-mode " 
        "0 -outdir results/{wildcards.library}/qualimap_{wildcards.library}_{wildcards.chromosome} -outformat HTML"
