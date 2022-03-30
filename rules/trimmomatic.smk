rule trimmomatic:
    input:
        expand("{libdir}{library}_{type}_{strand}.fastq", libdir=config["libdir"], type=["mutated", "base"], strand=["R1", "R2"], allow_missing=True)
    output:
        temp(expand("results/{library}/{library}_{type}_{strand}.trimmed.paired.fastq", type=["mutated", "base"], strand=["R1", "R2"], allow_missing=True)),
        temp(expand("results/{library}/{library}_{type}_{strand}.trimmed.unpaired.fastq", type=["mutated", "base"], strand=["R1", "R2"], allow_missing=True)),
        temp(expand("results/{library}/trim_{type}_out.log", type=["mutated", "base"], allow_missing=True))
    log:
        "results/logs/trimmomatic_{library}.log"
    benchmark:
        "results/benchmarks/trimmomatic_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/trimmomatic/trimmomatic:latest"
    threads: 4
    resources:
        mem_mb=16000
    shell:
        "trimmomatic PE -threads {threads} {input[0]} {input[1]} {output[0]} {output[4]} {output[1]} {output[5]} "
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:15 2> {output[8]} ; "
        "trimmomatic PE -threads {threads} {input[2]} {input[3]} {output[2]} {output[6]} {output[3]} {output[7]} "
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:15 2> {output[9]}"
