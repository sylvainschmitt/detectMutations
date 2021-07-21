## Sylvain SCHMITT
## 21/07/2021

import pandas as pd
configfile: "config/config.dag.yml"

chromosomes_table = pd.read_table(config["refdir"] + "/" + config["reference"] + ".fa.fai",
                                  header = None, names = ["chr", "X2", "X3", "X4", "X5"])
chromosomes = list(chromosomes_table.chr)
lambda wildcards: chromosomes
# print(expand("{chromosome}", chromosome=chromosomes))

libraries = expand("C{repetition}", repetition=config["repetitions"])
libraries.extend(expand("B{branch}_T{tip}_L{repetition}", branch=config["branches"], tip=config["tips"], repetition=config["repetitions"]))
lambda wildcards: libraries
# print(expand("{lib}", lib=libraries))

rule all:
    input:
        expand("results/reference/{reference}_{chromosome}.fa", reference=config["reference"], chromosome=chromosomes), # ref
        expand("results/reads/{library}_R{strand}.trimmed.paired.fq", library=libraries, strand=["1", "2"]), # reads
        expand("results/alns/{library}_on_{chromosome}.md.cram", library=libraries, chromosome=chromosomes), # alns
        expand("results/mutations/B{branch}_T{tip}_on_{chromosome}.tip.{ext}", 
                branch=config["branches"], tip=config["tips"], chromosome=chromosomes, ext=["vcf", "tsv"]) # muts

# Rules #

## Reference & reads ##
include: "rules/samtools_faidx_split.smk"
include: "rules/bwa_index.smk"
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
include: "rules/bedtools_subtract.smk"
include: "rules/bedtools_intersect.smk"
include: "rules/strelka2tsv.smk"
