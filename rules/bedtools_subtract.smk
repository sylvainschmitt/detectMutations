rule bedtools_subtract:
    input:
        "results/mutations/B{branch}_T{tip}_L{repetition}_on_{chromosome}.raw.vcf",
        expand("{refdir}/{hz}", refdir=config["refdir"], hz=config["hz"])
    output:
        temp("results/mutations/B{branch}_T{tip}_L{repetition}_on_{chromosome}.nonhz.vcf")
    log:
        "results/logs/bedtools_subtract_B{branch}_T{tip}_L{repetition}_on_{chromosome}.log"
    benchmark:
        "results/benchmarks/bedtools_subtract_B{branch}_T{tip}_L{repetition}_on_{chromosome}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools subtract -header -a {input[0]} -b {input[1]} > {output}"