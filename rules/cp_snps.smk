rule cp_snps:
    input:
        expand("{refdir}{snps}", refdir=config["refdir"], snps=config["snps"])
    output:
        expand("results/reference/{snps}", snps=config["snps"])
    log:
        "results/logs/cp_snps.log"
    benchmark:
        "results/benchmarks/cp_snps.benchmark.txt"
    shell:
        "cp {input} {output}"
        