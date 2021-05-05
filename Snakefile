## Sylvain SCHMITT
## 28/04/2021

configfile: "config/config.yml"

rule all:
    input:
        ## index ##
        # multiext("results/raw_data/reference/reference", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".dict"),
        ## reads qc ##
        "results/qc/reads/multiqc_report.html", 
        ## mapping ##
        expand("results/alignments/bwa/{sample}{ext}", sample=[config["base"], config["mutated"]], ext=[".md.bam", ".md.bam.bai", ".bam.metrics", ".md.bam.stats.out", "/qualimapReport.html"]),
        ## germline variant calling ##
        expand("results/germline/{caller}/{sample}.vcf", caller=["gatk", "freebayes"], sample=[config["base"]]),
        ## somatic variant calling ##
        expand("results/somatic/{caller}/{mutated}_vs_{base}.vcf", caller=["freebayes", "gatk"], base=[config["base"]], mutated=[config["mutated"]]),
        ## calls qc ##
        expand("results/somatic/{caller}/{mutated}_vs_{base}{ext}", caller=["freebayes", "gatk"], base=[config["base"]], 
                mutated=[config["mutated"]], ext=[".bcftools.stats.out", ".vcf.FILTER.summary", ".vcf.TsTv.count", ".vcf.TsTv.qual"])
        
# Rules #

# inspired from https://github.com/nf-core/sarek/blob/master/main.nf

## Index ##

include: "rules/mv_data.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/gatk_dict.smk"

## Reads QC ##
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"

## Trimming ##

# to be added

## Mapping ##

include: "rules/bwa_mem.smk"
include: "rules/samtools_sort.smk"
include: "rules/samtools_index.smk"
include: "rules/gatk_markduplicates.smk"

## Recalibration ##

# No know site, skipping BQSR for the moment
# include: "rules/gatk_baserecalibrator.smk"
# include: "rules/gatk_applybqsr.smk"

## Alignments QC ##

include: "rules/samtools_stats.smk"
include: "rules/qualimap.smk"

## Germline variant calling ##

# gatk haplotypecaller + genotypegvcfs, sentieon DNAseq + DNAscope, strelka, manta, tiddit, freebayes 
# first with gatk and freebayes

### freebayes ###

include: "rules/freebayes_germline.smk"

### GATK  ###

include: "rules/gatk_haplotypecaller.smk"
include: "rules/gatk_genotypegvcfs.smk"

## Somatic variant calling ##

# gatk mutect2, freebayes, strelka2, manta, ascat, control-freec, msisensor
# focus first on gatk and freebayes

### freebayes ###

include: "rules/freebayes_somatic.smk"

### GATK Mutect2 ###

include: "rules/gatk_mutect2.smk"
# include: "rules/gatk_mutect2.smk"

## Calls QC ##

include: "rules/bcftools_stats.smk"
include: "rules/vcftools_stats.smk"


## Annotation ##