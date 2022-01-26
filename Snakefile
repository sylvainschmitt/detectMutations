configfile: "config/config.yml"

rule all:
    input:
        expand("results/mutations/{comp}.raw.sql", comp=config["comps"]), # mut leaf
        "results/multiqc_report.html" #qc

# Rules #

## Reference ##
include: "rules/cp_reference.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"

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
include: "rules/mosdepth_regions.smk"

## Mutations ##
include: "rules/strelka2.smk"
include: "rules/strelka2tsv.smk"
include: "rules/strelka2sql.smk"

## QC ##
include: "rules/multiqc.smk"
