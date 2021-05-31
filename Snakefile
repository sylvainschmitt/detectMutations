## Sylvain SCHMITT
## 28/04/2021

configfile: "config/config.experiment.yml"
libraries, = glob_wildcards(config["libdir"] + "/{library}_mutated_R1.fastq")

rule all:
    input:
        ## mutations ##
        expand("results/mutations/{library}_{caller}.vcf", library=libraries, 
                caller=["mutect2", "freebayes", "gatk", "strelka2", "manta", "varscan", "somaticsniper", "muse"]),
        ## qc ##
        expand("results/{library}/multiqc_report.html", library=libraries)

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
include: "rules/trimmomatic.smk"

## Alignments ##
include: "rules/bwa_mem.smk"
include: "rules/samtools_sort.smk"
include: "rules/samtools_index.smk"
include: "rules/gatk_markduplicates.smk"
include: "rules/samtools_index_md.smk"
include: "rules/samtools_stats.smk"
include: "rules/samtools_mpileup.smk"
include: "rules/qualimap.smk"

## Mutations ##

### GATK Mutect2 ###
include: "rules/gatk_mutect2.smk"

### freebayes ###
include: "rules/freebayes.smk"
include: "rules/bedtools_substract.smk"

### GATK HaplotypeCaller ###
include: "rules/gatk_haplotypecaller.smk"
include: "rules/gatk_genotypegvcfs.smk"
include: "rules/bedtools_substract.smk"

## Strelka2
include: "rules/strelka2.smk"

## Manta
include: "rules/manta.smk"

## VarScan2
# not working
# include: "rules/varscan2.smk"

## VarScan
include: "rules/varscan.smk"
include: "rules/varscan2vcf.smk"

## Somatic Sniper
include: "rules/somaticsniper.smk"
include: "rules/bedtools_substract.smk"

## CaVEMan
include: "rules/caveman.smk"

## MuSe
include: "rules/muse.smk"

## RADIA
include: "rules/radia.smk"

## results ##
include: "rules/cp_vcfs.smk"
include: "rules/multiqc.smk"
