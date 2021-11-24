## Sylvain SCHMITT
## 21/07/2021

import pandas as pd
configfile: "config/config.dag.yml"

libraries, = glob_wildcards(config["libdir"] + "/{library}_1.fq.gz")

rule all:
    input:
        ## reads ##
        expand("results/{library}/{library}_{strand}.fastqc.{ext}", library=libraries,
                strand=["1", "2"], ext=["html", "zip"], allow_missing=True)
        # expand("results/reference/{reference}.fa", reference=config["reference"]), # ref
        # expand("results/reads/{library}_R{strand}.trimmed.paired.fq", library=libraries, strand=["1", "2"]), # reads
        # expand("results/alns/{library}.md.cram", library=libraries), # alns
        # expand("results/mutations/B{branch}_T{tip}.tip.{ext}", branch=config["branches"], tip=config["tips"], ext=["vcf", "tsv"]) # muts

# Rules #

## Reference ##
# include: "rules/cp_reference.smk"
# include: "rules/bwa_index.smk"

## Reads ##
include: "rules/trimmomatic.smk"
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"

## Alignments ##
# include: "rules/bwa_mem.smk"
# include: "rules/samtools_view.smk"
# include: "rules/samtools_sort.smk"
# include: "rules/samtools_index.smk"
# include: "rules/gatk_markduplicates.smk"
# include: "rules/samtools_view_md.smk"
# include: "rules/samtools_index_md.smk"

## Mutations ##
# include: "rules/strelka2.smk"
# include: "rules/bedtools_subtract.smk"
# include: "rules/bedtools_intersect.smk"
# include: "rules/strelka2tsv.smk"
