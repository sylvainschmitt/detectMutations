## Sylvain SCHMITT
## 28/04/2021

configfile: "config/config.yml"

rule all:
    input:
        # multiext("results/raw_data/reference/reference", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict"),
        "results/qc/reads/multiqc_report.html",
        expand("results/alignments/bwa/{sample}{ext}", sample=[config["base"], config["mutated"]], ext=[".md.bam", ".md.bam.bai", ".bam.metrics", ".md.bam.stats.out", "/qualimapReport.html"]),
        expand("results/germline/{caller}/{sample}.vcf", caller=["gatk", "freebayes"], sample=[config["base"]])
        
# Rules

# inspired from https://github.com/nf-core/sarek/blob/master/main.nf

## index
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/gatk_dict.smk"
## read QC
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
## trimming
# to be added
## mapping
include: "rules/bwa_mem.smk"
include: "rules/samtools_sort.smk"
include: "rules/samtools_index.smk"
include: "rules/gatk_markduplicates.smk"
## recalibration
# No know site, skipping BQSR for the moment
# include: "rules/gatk_baserecalibrator.smk"
# include: "rules/gatk_applybqsr.smk"
## alignments QC
include: "rules/samtools_stats.smk"
include: "rules/qualimap.smk"
## germline variant calling
# gatk haplotypecaller + genotypegvcfs, sentieon DNAseq + DNAscope, strelka, manta, tiddit, freebayes 
# focus first on gatk (maybe free abyes in second)
include: "rules/gatk_haplotypecaller.smk"
include: "rules/gatk_genotypegvcfs.smk"
include: "rules/freebayes.smk"
## somatic variant calling
# gatk mutect2, freebayes, strelka2, manta, ascat, control-freec, msisensor
# focus first on mutect2
# include: "rules/gatk_mutect2.smk"
# qc with bcftools stats and vcftools --TsTv-by-count --TsTv-by-qual -FILTER-summary 
