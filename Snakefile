## Sylvain SCHMITT
## 26/07/2021

import pandas as pd

configfile: "config/config.bordeaux.yml"

libraries, = glob_wildcards(config["libdir"] + "/{library}_1.fastq.gz")

chromosomes_table = pd.read_table(config["refdir"] + "/" + config["reference"] + ".fa.fai",
                                  header = None, names = ["chr", "X2", "X3", "X4", "X5"])
chromosomes = list(chromosomes_table.chr)
lambda wildcards: chromosomes
# print(expand("{chromosome}", chromosome=chromosomes))

rule all:
    input:
        expand("results/{library}/{library}_{strand}.trimmed.paired.fastq.gz", 
                strand=["1", "2"], library=libraries, chromosome=chromosomes) # trim


# Rules #

## Reference ##
include: "rules/samtools_faidx_split.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/gatk_dict.smk"

## Reads ##
include: "rules/cp_reads.smk"
include: "rules/trimmomatic.smk"

## Alignments ##
include: "rules/bwa_mem.smk"
include: "rules/samtools_view.smk"
include: "rules/samtools_sort.smk"
include: "rules/samtools_index.smk"
include: "rules/gatk_markduplicates.smk"
include: "rules/samtools_view_md.smk"
include: "rules/samtools_index_md.smk"

## Mutations ##
include: "rules/strelka2.smk"
include: "rules/strelka2tsv.smk"
include: "rules/strelka2sql.smk"
include: "rules/gatk_haplotypecaller.smk"
include: "rules/gatk_genotypegvcfs.smk"
include: "rules/gatk2tsv.smk"
include: "rules/gatk2sql.smk"
