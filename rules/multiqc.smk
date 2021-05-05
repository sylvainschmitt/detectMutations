rule multiqc:
    input:
        expand("results/qc/reads/{sample}_{strand}_fastqc.{ext}", sample=[config["base"], config["mutated"]], strand=["R1","R2"], ext=["html", "zip"])
    output:
        "results/qc/reads/multiqc_report.html"
    log:
        "results/logs/multiqc.log"
    benchmark:
        "results/benchmarks/multiqc.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/multiqc/multiqc:latest"
    shell:
        "multiqc results/qc/reads/ -o results/qc/reads/"

# could include bam reports from samtools stats