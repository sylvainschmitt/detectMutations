rule muse:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}_snps.vcf", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}.fa{ext}", 
               refdir=config["refdir"], reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}_snps.vcf.idx", refdir=config["refdir"], reference=config["reference"], allow_missing=True)
    output:
        "results/{lib}_REP{REP,\d+}/muse/{lib}_REP{REP}.vcf"
    log:
        "results/logs/muse_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/muse_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "docker://opengenomics/muse"
    threads: 4
    resources:
        mem_mb=16000
    shell:
        "muse.py -f {input[0]} --tumor-bam {input[1]} --normal-bam {input[2]} -O {output} -n {threads} -w results/{wildcards.lib}_REP{wildcards.REP}/muse/"
        