rule bcftools_call:
    input:
        "results/{library}/{library}_mutated.mpileup"
    output:
        temp("results/{library}/bcftools/{library}.unfiltered.vcf")
    log:
        "results/logs/bcftools_call_{library}.log"
    benchmark:
        "results/benchmarks/bcftools_call_{library}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest"
    shell:
        "bcftools call {input} -mv -Ov -o {output}"
