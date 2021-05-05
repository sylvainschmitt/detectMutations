sample=[config["base"], config["mutated"]]

rule fgbio_fastqtobam:
    input:
        expand("results/raw_data/reads/{sample}_{strand}.fastq", strand=["R1","R2"], allow_missing=True),
        expand("results/raw_data/reads/{sample}_{strand}.str", strand=["R1","R2"], allow_missing=True)
    output:
        "${idSample}_umi_converted.bam" 
    log:
        "results/logs/fgbio_fastqtobam_{sample}.log"
    benchmark:
        "results/benchmarks/fgbio_fastqtobam_{sample}.benchmark.txt"
    singularity: 
        "docker://zlskidmore/fgbio"
    shell:
        "fgbio --tmp-dir=${PWD}/tmp FastqToBam -i {input[0]} {input[1]} -o {output} "
        "--read-structures {input[2]} {input[3]} --sample {wildcards.sample} --library {wildcards.sample}"
