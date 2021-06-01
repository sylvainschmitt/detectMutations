## Sylvain SCHMITT
## 28/04/2021

import pandas as pd

configfile: "config/config.swiss.yml"

libraries, = glob_wildcards(config["libdir"] + "/{library}_1.fastq.gz")

chromosomes_table = pd.read_table(config["refdir"] + "/" + config["reference"] + ".fa.fai",
                                  header = None, names = ["chr", "X2", "X3", "X4", "X5"])
chromosomes = list(chromosomes_table.chr)
lambda wildcards: chromosomes
# print(expand("{chromosome}", chromosome=chromosomes))

rule all:
    input:
        # expand("results/reference/{reference}_{chromosome}.{ext}", 
        #         reference=config["reference"], ext=["fa", "fa.fai", "dict", "fa.amb"], chromosome=chromosomes),
        ## alignments ##
        expand("results/{library}/{library}_{chromosome}_md.bam", library=libraries, chromosome=chromosomes),
        ## qc ##
        "results/multiqc_report.html"

# Rules #

## Reference ##
include: "rules/cp_reference.smk"
include: "rules/samtools_faidx_split.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/gatk_dict.smk"

## Reads ##
include: "rules/cp_reads.smk"
include: "rules/fastqc.smk"
include: "rules/trimmomatic.smk"

## Alignments ##
include: "rules/bwa_mem.smk"
include: "rules/samtools_sort.smk"
include: "rules/samtools_index.smk"
include: "rules/gatk_markduplicates.smk"
include: "rules/samtools_index_md.smk"
include: "rules/samtools_stats.smk"
include: "rules/qualimap.smk"

## Mutations ##
# include: "rules/manta.smk"
# include: "rules/strelka2.smk"

## qc ##
include: "rules/multiqc.smk"
