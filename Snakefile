## Sylvain SCHMITT
## 27/07/2021

configfile: "config/config.swiss.yml"

libraries, = glob_wildcards(config["libdir"] + "/{library}_1.fastq.gz")

rule all:
    input:
        # expand("results/{library}/{library}.md.cram", library=libraries), # aln
        expand("results/mutations/{vcfs}_{caller}.vcf", vcfs=config["vcfs"], caller=["strelka2", "gatk"]) # mut raw vcf
        # expand("results/mutations/{chromosome}_{caller}.tsv", vcfs=config["vcfs"], chromosome=chromosomes, caller=["gatk"]), # mut raw tsv
        # "results/napoleon_mutations.tsv",
        # expand("results/{caller}_raw.sql", caller=["strelka2", "gatk"])


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
# include: "rules/strelka2tsv.smk"
# include: "rules/strelka2sql.smk"
include: "rules/gatk_haplotypecaller.smk"
include: "rules/gatk_cnnscorevariants.smk"
include: "rules/gatk_filtervarianttranches.smk"
# include: "rules/gatk2tsv.smk"
# include: "rules/gatk2sql.smk"
