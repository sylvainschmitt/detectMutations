rule strelka2:
    input:
        expand("{refdir}{reference}_REP{REP}.fa", refdir=config["refdir"], reference=config["reference"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{lib}_REP{REP}/{lib}_REP{REP}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("{refdir}{reference}_REP{REP}.fa{ext}", 
               refdir=config["refdir"], reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True)
    output:
        temp("results/{lib}_REP{REP,\d+}/strelka2/{lib}_REP{REP}.unfiltered.vcf")
    log:
        "results/logs/strelka2_{lib}_REP{REP}.log"
    benchmark:
        "results/benchmarks/strelka2_{lib}_REP{REP}.benchmark.txt"
    singularity: 
        "docker://quay.io/wtsicgp/strelka2-manta"
    threads: 4
    resources:
        mem_mb=16000
    shell:
        "rm -r results/{wildcards.lib}_{wildcards.REP}/strelka2 ; "
        "configureStrelkaSomaticWorkflow.py "
        "--normalBam {input[2]} "
        "--tumorBam {input[1]} "
        "--referenceFasta {input[0]} "
        "--runDir results/{wildcards.lib}_{wildcards.REP}/strelka2/ ; "
        "results/{wildcards.lib}_{wildcards.REP}/strelka2/runWorkflow.py -m local -j {threads} ; "
        "zcat results/{wildcards.lib}_{wildcards.REP}/strelka2/results/variants/somatic.snvs.vcf.gz > {output}"
        
