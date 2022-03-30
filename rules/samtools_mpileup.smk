rule samtools_mpileup:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}.fa{ext}", refdir=config["refdir"], 
                reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        expand("results/{lib}_REP{REP,\d+}/{lib}_REP{REP}_{type}.mpileup", type=["mutated", "base"], allow_missing=True)
    log:
        "results/logs/samtools_index_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/samtools_index_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools mpileup -f {input[0]} {input[1]} > {output[0]} ; "
        "samtools mpileup -f {input[0]} {input[2]} > {output[1]}"
