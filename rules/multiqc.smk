rule multiqc:
    input:
        expand("results/{library}/{library}_{type}_{strand}.raw_fastqc.{ext}", type=["mutated", "base"], strand=["R1", "R2"], ext=["html", "zip"], allow_missing=True),
        expand("results/{library}/trim_{type}_out.log", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.md.bam.stats.out", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/qualimap_{type}/qualimapReport.html", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.bam.metrics", type=["mutated", "base"], allow_missing=True)
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
        "rm -r results/{wildcards.library}/qualimap_mutated ; "
        "rm -r results/{wildcards.library}/qualimap_base ; "

# Can incolude outputs from:
# VarScan2
# VCFTools TsTv-by-count , TsTv-by-qual, TsTv-summary 
