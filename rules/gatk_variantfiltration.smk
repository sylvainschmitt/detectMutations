rule gatk_variantfiltration:
    input:
        "results/mutations/{tumor}_vs_{normal}_gatk.snps.vcf"
    output:
        "results/mutations/{tumor}_vs_{normal}_gatk.filtered.vcf"
    log:
        "results/logs/gatk_variantfiltration_{tumor}_vs_{normal}.log"
    benchmark:
        "results/benchmarks/gatk_variantfiltration_{tumor}_vs_{normal}.benchmark.txt"
    singularity: 
        "docker://broadinstitute/gatk"
    threads: 1
    resources:
        mem_mb=10000
    params:
        max_mem = lambda wildcards, resources: resources.mem_mb
    shell:
        "gatk VariantFiltration --java-options \"-Xmx{params.max_mem}M -Xms1G -Djava.io.tmpdir=tmp\" -V {input} "
        "-filter \"QD < 2.0\" --filter-name \"QD2\" "
        "-filter \"QUAL < 30.0\" --filter-name \"QUAL30\" "
        "-filter \"SOR > 3.0\" --filter-name \"SOR3\" "
        "-filter \"FS > 60.0\" --filter-name \"FS60\" "
        "-filter \"MQ < 40.0\" --filter-name \"MQ40\" "
        "-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" "
        "-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" "
        "-O {output}"
