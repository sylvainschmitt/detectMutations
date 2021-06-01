rule multiqc:
    input:
        expand("results/{library}/{library}_{strand}.raw_fastqc.{ext}", library=libraries, strand=["1", "2"], ext=["html", "zip"]),
        expand("results/{library}/trim_out.log", library=libraries),
        expand("results/{library}/{library}_{chromosome}_md.bam.stats.out", chromosome=chromosomes, library=libraries),
        expand("results/{library}/qualimap_{library}_{chromosome}/qualimapReport.html", chromosome=chromosomes, library=libraries),
        expand("results/{library}/{library}_{chromosome}.bam.metrics", chromosome=chromosomes, library=libraries)
    output:
        "results/multiqc_report.html"
    log:
        "results/logs/multiqc.log"
    benchmark:
        "results/benchmarks/multiqc.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/multiqc/multiqc:latest"
    shell:
        "multiqc results/ -o results/ ; rm -r results/*/qualimap_* ; rm -r results/*/aln"
