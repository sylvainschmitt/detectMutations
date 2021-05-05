caller=["freebayes", "gatk"]
base=[config["base"]]
mutated=[config["mutated"]]

rule bcftools_stats:
    input:
        "results/somatic/{caller}/{mutated}_vs_{base}.vcf"
    output:
        "results/somatic/{caller}/{mutated}_vs_{base}.bcftools.stats.out"
    log:
        "results/logs/bcftools_stats_{caller}_{mutated}_vs_{base}.log"
    benchmark:
        "results/benchmarks/bcftools_stats_{caller}_{mutated}_vs_{base}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest"
    shell:
        "bcftools stats {input} > {output}"
