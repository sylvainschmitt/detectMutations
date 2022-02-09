rule bedtools_substract_hz:
    input:
        "results/mutations_{tissue}/{tumor}_vs_{base}.raw.vcf",
        "results/hz/shared_hz.vcf.gz"
    output:
        "results/mutations_{tissue}/{tumor}_vs_{base}.nonhz.vcf"
    log:
        "results/logs/bedtools_substract_hz__{tissue}_{tumor}_vs_{base}.log"
    benchmark:
        "results/benchmarks/bedtools_substract_hz__{tissue}_{tumor}_vs_{base}.benchmark.txt"
    singularity: 
        "oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest"
    shell:
        "bedtools subtract -header -a {input[0]} -b {input[1]} > {output}"
        