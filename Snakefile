configfile: "config/config.yml"
# configfile: "config/config.dag.yml"

intervals, = glob_wildcards(config["intervals"] + "/{interval}")

rule all:
    input:
        "results/hz/shared_hz.vcf.gz", # hz
        expand("results/mutations_cambium/{comps}.nonhz.vcf", comps=config["cambium_comp"]), # mut cambium
        expand("results/mutations_leaf/{tumor}_vs_{base}.nonhz.vcf", tumor=config["leaf"], base=config["cambium_ref"]), # mut leaf
        "results/multiqc_report.html" #qc

# Rules #

## Reference ##
include: "rules/cp_reference.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/gatk_dict.smk"
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
include: "rules/gatk_haplotypecaller.smk"
include: "rules/gatk_gathergvcfs.smk"
include: "rules/gatk_genomicsdbimport.smk"
include: "rules/gatk_genotypegvcfs.smk"
include: "rules/gatk_gathervcfs.smk"
include: "rules/bcftools_biallelic.smk"
include: "rules/gatk_snps.smk"
include: "rules/plink_nonmissing.smk"
include: "rules/bcftools_shared.smk"

## Mutations ##
include: "rules/strelka2.smk"
include: "rules/bedtools_subtract_hz.smk"
# include: "rules/strelka2tsv.smk"
# include: "rules/strelka2sql.smk"

## QC ##
include: "rules/multiqc.smk"
