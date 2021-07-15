rule blat:
    input:
        expand("{refdir}/{reference}.fa", refdir=config["refdir"], reference=config["reference"]),
        "results/napoleon/napoleon_mutations.fa"
    output:
        temp("results/napoleon/napoleon_mutations.psl")
    log:
        "results/logs/blat.log"
    benchmark:
        "results/benchmarks/blat.benchmark.txt"
    shell:
        "~/Tools/blatSrc/bin/blat {input[0]} {input[1]} {output}"
