rule multiqc:
    input:
        expand("results/{library}/{library}_{strand}.raw_fastqc.{ext}", strand=["1", "2"], ext=["html", "zip"], allow_missing=True),
        "results/{library}/trim_out.log",
        "results/{library}/{library}.md.bam.stats.out",
        "results/{library}/qualimap/qualimapReport.html",
        "results/{library}/{library}.bam.metrics"
    output:
        "results/{library}/multiqc_report.html"
    log:
        "results/logs/multiqc_{library}.log"
    benchmark:
        "results/benchmarks/multiqc_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/multiqc/multiqc:latest"
    shell:
        "multiqc results/{wildcards.library}/ -o results/{wildcards.library} ; "
        "rm -r results/{wildcards.library}/qualimap"
