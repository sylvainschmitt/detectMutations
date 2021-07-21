rule trimmomatic:
    input:
        expand("{reads}/{library}_R{strand}.fq", reads=config["reads"] , strand=["1", "2"], allow_missing=True)
    output:
        expand("results/reads/{library}_R{strand}.trimmed.paired.fq", strand=["1", "2"], allow_missing=True),
        temp(expand("results/reads/{library}_R{strand}.trimmed.unpaired.fq", strand=["1", "2"], allow_missing=True)),
        temp("results/reads/{library}_trim_out.log")
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
        "trimmomatic PE -threads {threads} {input[0]} {input[1]} {output[0]} {output[2]} {output[1]} {output[3]} "
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:15 2> {output[4]}"