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
    	expand("results/{library}/{library}_{strand}.trimmed.paired.fastq.gz", library=libraries, strand=["1", "2"]), # reads
        expand("results/{library}/{library}_{chromosome}.md.cram", library=libraries, chromosome=chromosomes), # aln
        "results/multiqc_report.html", # qc
        expand("results/mutations/{vcfs}_on_{chromosome}_strelka2.vcf", vcfs=config["vcfs"], chromosome=chromosomes) # mut

# Rules #

## Reference ##
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
include: "rules/samtools_view.smk"
include: "rules/samtools_sort.smk"
include: "rules/samtools_index.smk"
include: "rules/gatk_markduplicates.smk"
include: "rules/samtools_view_md.smk"
include: "rules/samtools_index_md.smk"
include: "rules/samtools_stats.smk"
include: "rules/qualimap.smk"

## Mutations ##
include: "rules/strelka2.smk"

## QC ##
include: "rules/multiqc.smk"
