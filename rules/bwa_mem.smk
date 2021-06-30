rule bwa_mem:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}_{strand}.trimmed.paired.fastq", type=["mutated", "base"], strand=["R1", "R2"], allow_missing=True),
        expand("results/reference/{reference}{ext}", reference=config["reference"], ext=[".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict"])
    output:
        temp(expand("results/{library}/{library}_{type}.sam", type=["mutated", "base"], allow_missing=True))
    params:
        rg_mutated=r"@RG\tID:{library}_mutated\tSM:{library}_mutated",
        rg_base=r"@RG\tID:{library}_base\tSM:{library}_base"
    log:
        "results/logs/bwa_mem_{library}.log"
    benchmark:
        "results/benchmarks/bwa_mem_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    threads: 8
    resources:
        mem_mb=32000
    shell:
        "bwa mem -M -R '{params.rg_mutated}' -t {threads} {input[0]} {input[1]} {input[2]} > {output[0]} ; "
        "bwa mem -M -R '{params.rg_base}' -t {threads} {input[0]} {input[3]} {input[4]} > {output[1]} ; "
        