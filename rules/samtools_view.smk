rule samtools_view:
    input:
        "results/alns/{library}_on_{chromosome}.sam",
        expand("results/reference/{reference}_{chromosome}.fa", reference=config["reference"], allow_missing = True)
    output:
        temp("results/alns/{library}_on_{chromosome}.raw.cram")
    log:
        "results/logs/samtools_view_{library}_{chromosome}.log"
    benchmark:
        "results/benchmarks/samtools_view_{library}_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools view -C -T {input[1]} -f 1 -F 12 {input[0]} > {output}"
