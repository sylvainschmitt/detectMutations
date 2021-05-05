aligner=["bwa"]

rule mutect2:
    input:
        "results/raw_data/reference/reference.fa",
        expand("results/alignments/{aligner}/{sample}.sorted.bam", sample=config["mutated"], allow_missing=True),
        expand("results/alignments/{aligner}/{sample}.sorted.bam", sample=config["base"], allow_missing=True),
        "results/raw_data/reference/reference.dict"
    output:
        multiext("results/mutations/mutect2/{aligner}/mutations", ".vcf.gz", ".vcf.gz.stats")
    log:
        "results/logs/mutect2_{aligner}.log"
    benchmark:
        "results/benchmarks/mutect2_{aligner}.benchmark.txt"
    singularity: 
        "docker://alexcoppe/gatk"
    shell:
        "/.singularity.d/runscript -T MuTect2   -R {input[0]} -I:tumor {input[1]} -I:normal {input[2]} -o {output[0]}"
