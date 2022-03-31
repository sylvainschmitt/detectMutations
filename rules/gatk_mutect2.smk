rule gatk_mutect2:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}_snps.vcf", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}.{ext}", 
               refdir=config["refdir"], reference=config["reference"], ext="dict", allow_missing=True),
        expand("{refdir}{reference}_REP{REP}_snps.vcf.idx", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    output:
        temp("results/{lib}_REP{REP,\d+}/mutect2/{lib}_REP{REP}.unfiltered.vcf")
    log:
        "results/logs/gatk_mutect2_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/gatk_mutect2_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk Mutect2 -R {input[0]} -I {input[1]} -tumor {wildcards.lib}_REP{wildcards.REP}_mutated  -I {input[2]} -normal {wildcards.lib}_REP{wildcards.REP}_base "
        " --panel-of-normals {input[3]} -dont-use-soft-clipped-bases true -O {output}"
        