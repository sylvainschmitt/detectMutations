rule bwa_mem:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}_{strand}.trimmed.paired.fastq", type=["mutated", "base"], strand=["R1", "R2"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], 
                reference=config["reference"], ext=[".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict"], allow_missing=True)
    output:
        temp(expand("results/{lib}_REP{REP,\d+}/{lib}_REP{REP}_{type}.sam", type=["mutated", "base"], allow_missing=True))
    params:
        rg_mutated=r"@RG\tID:{lib}_REP{REP}_mutated\tSM:{lib}_REP{REP}_mutated",
        rg_base=r"@RG\tID:{lib}_REP{REP}_base\tSM:{lib}_REP{REP}_base"
    log:
        "results/logs/bwa_mem_{lib}_REP{REP,\d+}.log"
    benchmark:
        "results/benchmarks/bwa_mem_{lib}_REP{REP,\d+}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    threads: 8
    resources:
        mem_mb=32000
    shell:
        "bwa mem -M -R '{params.rg_mutated}' -t {threads} {input[0]} {input[1]} {input[2]} > {output[0]} ; "
        "bwa mem -M -R '{params.rg_base}' -t {threads} {input[0]} {input[3]} {input[4]} > {output[1]} ; "
        