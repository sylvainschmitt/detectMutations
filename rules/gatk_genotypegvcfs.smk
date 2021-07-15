rule gatk_genotypegvcfs:
    input:
        expand("results/reference/{reference}_{chromosome}.fa", 
                reference=config["reference"], allow_missing=True),
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_gatk.g.vcf",
        expand("results/reference/{reference}_{chromosome}.fa{ext}", 
                reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], allow_missing=True),
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_gatk.g.vcf.idx"
    output:
        "results/mutations/{tumor}_vs_{normal}_on_{chromosome}_gatk.vcf"
    log:
        "results/logs/gatk_genotypegvcfs_{tumor}_vs_{normal}_{chromosome}.log"
    benchmark:
        "results/benchmarks/gatk_genotypegvcfs_{tumor}_vs_{normal}_{chromosome}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 20
    resources:
        mem_mb=100000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk GenotypeGVCFs --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -R {input[0]} -V {input[1]} -O {output}"
        
        