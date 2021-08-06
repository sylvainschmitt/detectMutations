## Sylvain SCHMITT
## 26/07/2021

import pandas as pd

configfile: "config/config.bordeaux.yml"

libraries, = glob_wildcards(config["libdir"] + "/{library}_1.fastq.gz")
intervals, = glob_wildcards(config["refdir"] + "/intervals/{intervals}")

rule all:
    input:
        expand("results/{library}/{library}.md.cram", library=libraries),
        expand("results/{caller}_raw.sql", caller=["strelka2", "mutect2"])

# Rules #

## Reference & Reads ##
include: "rules/cp_reference.smk"
include: "rules/bwa_index.smk"
include: "rules/samtools_faidx.smk"
include: "rules/gatk_dict.smk"
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
include: "rules/strelka2tsv.smk"
include: "rules/strelka2sql.smk"
include: "rules/gatk_mutect2.smk"
include: "rules/gatk_gathervcfs.smk"
include: "rules/mutect2tsv.smk"
include: "rules/mutect2sql.smk"
