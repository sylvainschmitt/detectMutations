caller2mv=["mutect2", "freebayes", "gatk", "strelka2", "manta", "varscan", "somaticsniper", "muse", "octopus"]

rule cp_vcfs:
    input:
        "results/{library}/{caller2mv}/{library}.vcf"
    output:
        "results/mutations/{library}_{caller2mv}.vcf"
    log:
        "results/logs/cp_{library}_{caller2mv}.log"
    benchmark:
        "results/benchmarks/cp_{library}_{caller2mv}.benchmark.txt"
    shell:
        "cp {input} {output}"
        