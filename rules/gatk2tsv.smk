rule gatk2tsv:
    input:
        "results/mutations/{vcfs}_gatk.vcf"
    output:
        "results/mutations/{vcfs}_gatk.tsv"
    log:
        "results/logs/gatk2tsv_{vcfs}.log"
    benchmark:
        "results/benchmarks/gatk2tsv_{vcfs}.benchmark.txt"
    shell:
        "bcftools query -f '%CHROM %POS %AC %AF %BaseQRankSum %INFO/DP %FS %MQ %MQRankSum %QD %ReadPosRankSum %SOR "
        "[%AD %DP %GQ]\\n' "
        "{input} > {output}"
        