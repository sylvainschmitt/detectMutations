## Sylvain SCHMITT
## 28/04/2021

configfile: "config/config.swiss.yml"
libraries, = glob_wildcards(config["libdir"] + "/{library}_1.fastq.gz")

rule all:
    input:
        ## alignments ##
        expand("results/{library}/{library}.md.bam", library=libraries),
        ## qc ##
        expand("results/{library}/multiqc_report.html", library=libraries)

# Rules #

## Reference ##
include: "rules/cp_reference.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/gatk_dict.smk"

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
include: "rules/qualimap.smk"

## Mutations ##
# include: "rules/manta.smk"
# include: "rules/strelka2.smk"

## qc ##
include: "rules/multiqc.smk"
