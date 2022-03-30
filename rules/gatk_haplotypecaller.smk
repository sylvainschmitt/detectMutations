rule gatk_haplotypecaller:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        "results/{lib}_REP{REP}/{lib}_REP{REP}_mutated.md.bam",
        "results/{lib}_REP{REP}/{lib}_REP{REP}_mutated.md.bam.bai",
        expand("{refdir}{reference}_REP{REP}.fa{ext}", 
               refdir=config["refdir"], reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        "results/{lib}_REP{REP,\d+}/gatk/{lib}_REP{REP}.g.vcf"
    log:
        "results/logs/gatk_haplotypecaller_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/gatk_haplotypecaller_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    shell:
        "gatk HaplotypeCaller -R {input[0]} -I {input[1]} -O {output} -ERC GVCF"
