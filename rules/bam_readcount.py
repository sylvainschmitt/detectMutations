rule bam_readcount:
    input:
        "results/alns/{library}_on_{reference}.md.cram",
        "results/reference/{reference}.fa",
        "data/angela.bed"
    output:
        "results/count/{library}_on_{reference}.tsv"
    log:
        "results/logs/bam_readcount_{library}_on_{reference}.log"
    benchmark:
        "results/benchmarks/bam_readcount_{library}_on_{reference}.benchmark.txt"
    singularity: 
        "docker://mgibio/bam-readcount"
    shell:
        "bam-readcount -f {input[1]} -l {input[2]} {input[0]} > {output}"
