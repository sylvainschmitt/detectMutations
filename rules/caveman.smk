rule caveman:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        "results/{library}/caveman/{library}.vcf"
    log:
        "results/logs/caveman_{library}.log"
    benchmark:
        "results/benchmarks/caveman_{library}.benchmark.txt"
    singularity: 
        "docker://leukgen/docker-caveman:v1.0.0"
    shell:
        "awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS \"0\"}}' {input[0]}.fai > results/{wildcards.library}/caveman/ignore.bed ; "
        "/.run caveman -o results/{wildcards.library}/caveman -r {input[0]}.fai -tb {input[1]} -nb {input[2]} "
        "-ig results/{wildcards.library}/caveman/ignore.bed -s {wildcards.library} -sa {wildcards.library} -st genome"
        