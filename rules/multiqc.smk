rule multiqc:
    input:
        expand("results/{library}/{library}_{strand}_fastqc.{ext}", library=config["samples"],
                strand=["1", "2"], ext=["html", "zip"]),
        expand("results/{library}/trim_out.log", library=config["samples"]),
        expand("results/alns/{library}_on_{reference}.md.cram.stats", library=config["samples"], reference=config["references"]),
        expand("results/alns/{library}_on_{reference}.mosdepth.global.dist.txt", library=config["samples"], reference=config["references"])
    output:
        "results/multiqc_report.html"
    log:
        "results/logs/multiqc.log"
    benchmark:
        "results/benchmarks/multiqc.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/multiqc/multiqc:latest"
    shell:
        "multiqc "
        "results/*/*_fastqc.zip "
        "results/*/trim_out.log "
        "results/alns/*.md.cram.stats "
        "results/alns/*.mosdepth.global.dist.txt "
        "-o results/"
        