## Sylvain SCHMITT
## 28/04/2021

configfile: "config/config.yml"
libraries, = glob_wildcards(config["libdir"] + "/{library}_R1.fastq")

rule all:
    input:
        ## reads ##
        expand("results/{library}/qc/multiqc_report.html", library=libraries),
        ## alignments ##
        expand("results/{library}/qc/qualimap/qualimapReport.html", library=libraries),
        expand("results/{library}/{library}.md.{ext}", library=libraries, ext=["bam", "bam.bai"]),
        ## mutations ##
        expand("results/{library}/{caller}/{library}.vcf", library=libraries, caller=["mutect2"])
        
        # expand("results/alignments/bwa/{sample}{ext}", sample=[config["base"], config["mutated"]], ext=[".md.bam", ".md.bam.bai", ".bam.metrics", ".md.bam.stats.out", "/qualimapReport.html"]),
        ## germline variant calling ##
        # expand("results/germline/{caller}/{sample}.vcf", caller=["gatk", "freebayes"], sample=[config["base"]]),
        ## somatic variant calling ##
        # expand("results/somatic/{caller}/{mutated}_vs_{base}.vcf", caller=["freebayes", "gatk"], base=[config["base"]], mutated=[config["mutated"]])


# Rules #

## Reference ##
include: "rules/cp_reference.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/gatk_dict.smk"
include: "rules/cp_snps.smk"
include: "rules/gatk_idx.smk"

## Reads ##
include: "rules/cp_reads.smk"
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
include: "rules/trimmomatic.smk"

## Alignments ##
include: "rules/bwa_mem.smk"
include: "rules/samtools_sort.smk"
include: "rules/samtools_index.smk"
include: "rules/gatk_markduplicates.smk"
include: "rules/samtools_index_md.smk"
include: "rules/qualimap.smk"

## Mutations ##

### GATK Mutect2 ###
include: "rules/gatk_mutect2.smk"

### GATK  ###
# include: "rules/gatk_haplotypecaller.smk"
# include: "rules/gatk_genotypegvcfs.smk"

### freebayes ###
# include: "rules/freebayes_somatic.smk"
