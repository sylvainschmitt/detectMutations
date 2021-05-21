rule varscan2:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        "results/{library}/varscan2/{library}.vcf"
    log:
        "results/logs/varscan2_{library}.log"
    benchmark:
        "results/benchmarks/varscan2_{library}.benchmark.txt"
    singularity: 
        "docker://jeltje/varscan2"
    shell:
        "awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS 5000}}' {input[0]}.fai > results/{wildcards.library}/varscan2/exome.bed ; "
        "awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS $2}}' {input[0]}.fai > results/{wildcards.library}/varscan2/exome.bed ; "
        "awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS \"0\"}}' {input[0]}.fai > results/{wildcards.library}/varscan2/centromere.bed ; "
        "/.run -c {input[2]} -t  {input[1]} -q {wildcards.library} -i {input[0]} "
        "-b results/{wildcards.library}/varscan2/centromere.bed -w results/{wildcards.library}/varscan2/exome.bed "
        "-s results/{wildcards.library}/varscan2 > {output}"
