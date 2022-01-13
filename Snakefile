configfile: "config/config.yml"
# configfile: "config/config.dag.yml"

rule all:
    input:
        expand("results/alns/{library}.md.cram", library=config["leaf"]), # alns
        # "results/leaf_nontrunk_mutations.sql", # mut
        # "results/trunk_raw_mutations.sql", # mut
        "results/multiqc_report.html" #qc

# Rules #

## Reference ##
include: "rules/cp_reference.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/bedtools_makewindows.smk"
include: "rules/bedtools_nuc.smk"
include: "rules/busco.smk"

## Reads ##
include: "rules/trimmomatic.smk"
include: "rules/fastqc.smk"

## Alignments ##
include: "rules/bwa_mem.smk"
include: "rules/samtools_view.smk"
include: "rules/samtools_sort.smk"
include: "rules/samtools_index.smk"
include: "rules/gatk_markduplicates.smk"
include: "rules/samtools_view_md.smk"
include: "rules/samtools_index_md.smk"
include: "rules/samtools_stats.smk"
include: "rules/mosdepth.smk"
include: "rules/mosdepth_regions.smk"

## Heterozygosity ##

## Mutations cambium ##
# include: "rules/strelka2.smk"
# include: "rules/bedtools_subtract.smk"
# include: "rules/strelka2tsv_leaf.smk"
# include: "rules/strelka2sql_leaf.smk"
# include: "rules/strelka2tsv_trunk.smk"
# include: "rules/strelka2sql_trunk.smk"

## Mutations leaf ##

## QC ##
include: "rules/multiqc.smk"
