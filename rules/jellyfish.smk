rule jellyfish:
    input:
        expand("results/{library}/{library}_{strand}.trimmed.paired.fq.gz", strand=["1", "2"], library=config["leaf_cov"])
    output:
        "results/hz/angela_21mer.jf",
        "results/hz/angela_21mer.histo"
    log:
        "results/logs/jellyfish.log"
    benchmark:
        "results/benchmarks/jellyfish.benchmark.txt"
    singularity: 
        "docker://quay.io/biocontainers/jellyfish:1.1.12--h6bb024c_1"
    threads: 4
    resources:
        mem_mb=40000
    params:
        libraries=lambda wildcards, input: "<(zcat " + ") <(zcat ".join(sorted(input)),
    shell:
        "jellyfish count -t {threads} -C -m 21 -s 8G -o {output[0]} {params.libraries}) ; "
        "jellyfish histo -t {threads} --high=1000000 -o {output[1]} {output[0]}"
        