rule jellyfish:
    input:
        expand("results/{library}/{library}_{strand}.trimmed.paired.fq.gz", strand=["1", "2"], library=config["hz"])
    output:
        "results/hz/angela_100mer.jf",
        "results/hz/angela_100mer.histo"
    log:
        "results/logs/jellyfish.log"
    benchmark:
        "results/benchmarks/jellyfish.benchmark.txt"
    singularity: 
        "docker://quay.io/biocontainers/jellyfish:1.1.12--h6bb024c_1"
    threads: 10
    resources:
        mem_mb=40000
    params:
        libraries=lambda wildcards, input: " < (zcat " + ") < (zcat ".join(sorted(input)),
        mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "jellyfish count -t {threads} -C -m 21 -s {params.mem}M -o {output[0]} {params.libraries}) ; "
        "jellyfish histo -t {threads} -o {output[1]} {output[0]}"
        