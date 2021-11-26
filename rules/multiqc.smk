rule multiqc:
    input:
        expand("results/{library}/{library}_{strand}_fastqc.{ext}", library=libraries,
                strand=["1", "2"], ext=["html", "zip"]),
        expand("results/{library}/trim_out.log", library=libraries)
    output:
        "results/multiqc_report.html"
    log:
        "results/logs/multiqc.log"
    benchmark:
        "results/benchmarks/multiqc.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/multiqc/multiqc:latest"
    shell:
        "multiqc results/*/*_fastqc.zip results/*/trim_out.log -o results/"