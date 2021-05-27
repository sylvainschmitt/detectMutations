rule trimmomatic:
    input:
        expand("results/{library}/{library}_{strand}.raw.fastq", strand=["1", "2"], allow_missing=True)
    output:
        temp(expand("results/{library}/{library}_{strand}.trimmed.paired.fastq", strand=["1", "2"], allow_missing=True)),
        temp(expand("results/{library}/{library}_{strand}.trimmed.unpaired.fastq", strand=["1", "2"], allow_missing=True)),
        temp("results/{library}/trim_out.log")
    log:
        "results/logs/trimmomatic_{library}.log"
    benchmark:
        "results/benchmarks/trimmomatic_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/trimmomatic/trimmomatic:latest"
    shell:
        "trimmomatic PE {input[0]} {input[1]} {output[0]} {output[2]} {output[1]} {output[3]} "
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:15 2> {output[4]}"