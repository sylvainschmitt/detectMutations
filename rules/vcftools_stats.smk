caller=["freebayes", "gatk"]
base=[config["base"]]
mutated=[config["mutated"]]

rule vcftools_stats:
    input:
        "results/somatic/{caller}/{mutated}_vs_{base}.vcf"
    output:
        temp("results/somatic/{caller}/{mutated}_vs_{base}.vcf.log"),
        "results/somatic/{caller}/{mutated}_vs_{base}.vcf.FILTER.summary",
        "results/somatic/{caller}/{mutated}_vs_{base}.vcf.TsTv.count",
        "results/somatic/{caller}/{mutated}_vs_{base}.vcf.TsTv.qual"
    log:
        "results/logs/vcftools_stats_{caller}_{mutated}_vs_{base}.log"
    benchmark:
        "results/benchmarks/vcftools_stats_{caller}_{mutated}_vs_{base}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/vcftools/vcftools:latest"
    shell:
        "vcftools --gzvcf {input} --TsTv-by-count --out {input} ; "
        "vcftools --gzvcf {input} --TsTv-by-qual --out {input} ; "
        "vcftools --gzvcf {input} --FILTER-summary --out {input}"
