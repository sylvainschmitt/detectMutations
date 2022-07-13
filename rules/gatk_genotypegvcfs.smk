rule gatk_genotypegvcfs:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        "results/{lib}_REP{REP}/gatk/{lib}_REP{REP,\d+}.g.vcf",
        expand("{refdir}{reference}_REP{REP}.fa{ext}", 
               refdir=config["refdir"], reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        "results/{lib}_REP{REP,\d+}/gatk/{lib}_REP{REP}.unfiltered.vcf"
    log:
        "results/logs/gatk_genotypegvcfs_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/gatk_genotypegvcfs_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk:4.2.6.1"
    shell:
        "gatk GenotypeGVCFs -R {input[0]} -V {input[1]} -O {output}"
        