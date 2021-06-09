rule octopus:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        temp("results/{library}/octopus/{library}.unfiltered.vcf")
    log:
        "results/logs/octopus_{library}.log"
    benchmark:
        "results/benchmarks/octopus_{library}.benchmark.txt"
    singularity: 
        "https://github.com/sylvainschmitt/singularity-octopus/releases/download/0.0.1/sylvainschmitt-singularity-octopus.latest.sif"
    shell:
        "octopus -R {input[0]} -I {input[2]}  {input[1]} -N {wildcards.library}_base -o {output}"