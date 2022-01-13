rule multiqc:
    input:
        directory(expand("results/reference/{reference}_busco", reference=config["reference"])),
        expand("results/reference/{reference}.gc", reference=config["reference"]),
        expand("results/{library}/{library}_{strand}_fastqc.{ext}", library=config["hz"],
                strand=["1", "2"], ext=["html", "zip"]),
        expand("results/{library}/trim_out.log", library=config["hz"]),
        expand("results/alns/{library}.md.cram.stats", library=config["hz"]),
        expand("results/alns/{library}.mosdepth.global.dist.txt", library=config["hz"]),
        expand("results/alns/{library}.mosdepth.region.dist.txt", library=config["leaf_cov"])
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
        "results/reference/*/short_summary.*.txt "
        "results/*/*_fastqc.zip "
        "results/*/trim_out.log "
        "results/alns/*.md.cram.stats "
        "results/alns/*.mosdepth.global.dist.txt "
        "results/alns/*.mosdepth.region.dist.txt "
        "-o results/"
        