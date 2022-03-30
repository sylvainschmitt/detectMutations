rule somaticsniper:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}.fa{ext}", 
               refdir=config["refdir"], reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        temp("results/{lib}_REP{REP,\d+}/somaticsniper/{lib}_REP{REP}.unfiltered.vcf")
    log:
        "results/logs/somaticsniper_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/somaticsniper_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "docker://lethalfang/somaticsniper:1.0.5.0"
    shell:
        "/opt/somatic-sniper/build/bin/bam-somaticsniper -q 25 -Q 15 -s 0.0001 -F vcf -f {input[0]} {input[1]}  {input[2]} {output}"