rule multiqc:
    input:
        expand("results/{library}/{library}_{strand}.raw_fastqc.{ext}", strand=["R1", "R2"], ext=["html", "zip"], allow_missing=True)
    output:
        "results/{library}/qc/multiqc_report.html"
    log:
        "results/logs/multiqc_{library}.log"
    benchmark:
        "results/benchmarks/multiqc_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/multiqc/multiqc:latest"
    shell:
        "multiqc results/{wildcards.library}/ -o results/{wildcards.library}/qc"

# could include bam reports from samtools stats
