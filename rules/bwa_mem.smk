sample=[config["base"], config["mutated"]]

rule bwa_mem:
    input:
        "results/raw_data/reference/reference.fa",
        "results/raw_data/reads/{sample}_R1.fastq",
        "results/raw_data/reads/{sample}_R2.fastq",
        multiext("results/raw_data/reference/reference", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict")
    output:
        temp("results/alignments/bwa/{sample}.sam")
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "results/logs/bwa_mem_{sample}.log"
    benchmark:
        "results/benchmarks/bwa_mem_{sample}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input[0]} {input[1]} {input[2]} > {output}"
        
# Help

# https://github.com/nf-core/sarek/blob/7ccfb36509d8380946c048065129b29d69c8443b/main.nf#L1177

# bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}

# bwa mem -K 100000000 -M -R "${readGroup}" -t {threads} ${input.fasta}
# samtools sort --threads {threads} {input.bam} -m 2G - > ${idSample}_${idRun}.bam
