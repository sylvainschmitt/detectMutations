rule radia:
    input:
        expand("results/reference/{reference}.fa", reference=config["reference"]),
        expand("results/{library}/{library}_{type}.md.bam", type=["mutated", "base"], allow_missing=True),
        expand("results/{library}/{library}_{type}.md.bam.bai", type=["mutated", "base"], allow_missing=True),
        expand("results/reference/{reference}.fa{ext}", reference=config["reference"], ext=[".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        "results/{library}/radia/{library}.vcf"
    log:
        "results/logs/radia_{library}.log"
    benchmark:
        "results/benchmarks/radia_{library}.benchmark.txt"
    singularity: 
        "docker://opengenomics/radia"
    shell:
        "/opt/radia.py --patientId {wildcards.library} -n {input[2]} -t {input[1]} -f {input[0]} -o {output}"
        
     