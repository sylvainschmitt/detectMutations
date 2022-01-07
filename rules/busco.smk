rule busco:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"])
    output:
        directory(expand("results/reference/{reference}_busco", reference=config["reference"]))
    log:
        "results/logs/busco.log"
    benchmark:
        "results/benchmarks/busco.benchmark.txt"
    singularity:
        "docker://ezlabgva/busco:v5.2.2_cv2"
    threads: 10
    resources:
        mem_mb=100000
    shell:
        "busco -i {input} -o {config[reference]}_busco -m genome -l viridiplantae --cpu {threads} --download_path tmp ; "
        "mv {config[reference]}_busco {output}"
        