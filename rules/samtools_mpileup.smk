rule samtools_mpileup:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        expand("results/{library}/{library}_{type}.mpileup", type=["mutated", "base"], allow_missing=True)
    log:
        "results/logs/samtools_index_{library}.log"
    benchmark:
        "results/benchmarks/samtools_index_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools mpileup -f {input[0]} {input[1]} > {output[0]} ; "
        "samtools mpileup -f {input[0]} {input[2]} > {output[1]}"
