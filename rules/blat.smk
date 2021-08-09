rule blat:
    input:
       expand("results/reference/{reference}.fa", reference=config["reference"]),
       "results/3P/mutations.fa"
    output:
        temp("results/3P/mutations.psl")
    log:
        "results/logs/blat.log"
    benchmark:
        "results/benchmarks/blat.benchmark.txt"
    shell:
        "~/Tools/blatSrc/bin/blat {input[0]} {input[1]} {output}"
        