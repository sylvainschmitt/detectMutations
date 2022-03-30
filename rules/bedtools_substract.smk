rule bedtools_substract:
    input:
        "results/{lib}_REP{REP}/{caller}/{lib}_REP{REP}.unfiltered.vcf",
        expand("{refdir}{reference}_REP{REP}_snps.vcf", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    output:
        "results/{lib}_REP{REP,\d+}/{caller}/{lib}_REP{REP}.vcf"
    log:
        "results/logs/bedtools_substract_{lib}_REP{REP}_{caller}.log"
    benchmark:
        "results/benchmarks/bedtools_substract_{lib}_REP{REP}_{caller}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools subtract -header -a {input[0]} -b {input[1]} > {output}"
