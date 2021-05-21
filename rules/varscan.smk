rule varscan:
    input:
        expand("results/{library}/{library}_{type}.mpileup", type=["mutated", "base"], allow_missing=True)
    output:
        temp("results/{library}/varscan/{library}.snp")
    log:
        "results/logs/varscan_{library}.log"
    benchmark:
        "results/benchmarks/varscan_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/varscan/varscan:latest"
    shell:
        "/.run somatic {input[1]} {input[0]} results/{wildcards.library}/varscan/{wildcards.library}"
