rule samtools_faidx:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    output:
        expand("{refdir}{reference}_REP{REP,\d+}.fa.fai", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    log:
        "results/logs/samtools_faidx_REP{REP}.log"
    benchmark:
        "results/benchmarks/samtools_faidx_REP{REP}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest"
    shell:
        "samtools faidx {input}"
