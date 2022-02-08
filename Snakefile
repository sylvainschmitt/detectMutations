configfile: "config/config.yml"

rule all:
    input:
        # expand("results/mutations/{comp}_on_{reference}.raw.vcf", comp=config["comps"], reference=config["references"]), # mut
        "results/mutations/mutations.robust.tsv",
        # "results/multiqc_report.html" #qc,
         "results/report.html"

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
include: "rules/circos_cov.smk"

## Mutations ##
include: "rules/strelka2.smk"
include: "rules/strelka2tsv.smk"
# include: "rules/strelka2sql.smk"
include: "rules/filter_mutations.smk"

## Report ##
include: "rules/multiqc.smk"
include: "rules/report.smk"
