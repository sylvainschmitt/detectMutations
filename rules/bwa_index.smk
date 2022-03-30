rule bwa_index:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    output:
        expand("{refdir}{reference}_REP{REP,\d+}.fa{ext}", 
                refdir=config["refdir"], reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    log:
        "results/logs/bwa_index_REP{REP}.log"
    benchmark:
        "results/benchmarks/bwa_index_REP{REP}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest"
    shell:
        "bwa index {input}"
