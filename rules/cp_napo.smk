rule cp_napo:
    input:
        expand("{libdir}/napoleon.fa", libdir=config["libdir"]),
        expand("{libdir}/napoleon_mutations.tsv", libdir=config["libdir"])
    output:
        temp("results/napoleon/napoleon.fa"),
        temp("results/napoleon/napoleon_mutations.tsv")
    log:
        "results/logs/cp_napo.log"
    benchmark:
        "results/benchmarks/cp_napo.benchmark.txt"
    shell:
        "cp {input[0]} {output[0]} ; cp {input[1]} {output[1]}"
        