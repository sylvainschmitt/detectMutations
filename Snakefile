## Sylvain SCHMITT
## 27/07/2021

configfile: "config/config.swiss.yml"

libraries, = glob_wildcards(config["libdir"] + "/{library}_1.fastq.gz")
intervals, = glob_wildcards(config["libdir"] + "/intervals/{intervals}")

rule all:
    input:
        "results/napoleon_mutations.tsv",
        expand("results/{caller}_raw.sql", caller=["strelka2", "gatk"])


# Rules #

## Reference ##
include: "rules/cp_reference.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/gatk_dict.smk"

## Napoleon ##
include: "rules/cp_napo.smk"
include: "rules/samtools_faidx_napo.smk"
include: "rules/napomutations2bed.smk"
include: "rules/bedtools_getfasta.smk"
include: "rules/blat.smk"
include: "rules/psl2pos.smk"

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
include: "rules/gatk_gathervcfs.smk"
include: "rules/gatk_selectsnps.smk"
include: "rules/gatk_variantfiltration.smk"
include: "rules/gatk_selectfiltered.smk"
include: "rules/gatk2tsv.smk"
include: "rules/gatk2sql.smk"
