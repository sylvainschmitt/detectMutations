rule gatk_idx:
    input:
        expand("{refdir}{reference}_REP{REP}_snps.vcf", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    output:
        expand("{refdir}{reference}_REP{REP,\d+}_snps.vcf.idx", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    log:
        "results/logs/gatk_idx_REP{REP}.log"
    benchmark:
        "results/benchmarks/gatk_idx_REP{REP}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk:4.2.6.1"
    shell:
        "gatk IndexFeatureFile -I {input}"
