rule somaticsniper:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        temp("results/{library}/somaticsniper/{library}.unfiltered.vcf")
    log:
        "results/logs/somaticsniper_{library}.log"
    benchmark:
        "results/benchmarks/somaticsniper_{library}.benchmark.txt"
    singularity: 
        "docker://lethalfang/somaticsniper:1.0.5.0"
    shell:
        "/opt/somatic-sniper/build/bin/bam-somaticsniper -q 25 -Q 15 -s 0.0001 -F vcf -f {input[0]} {input[1]}  {input[2]} {output}"