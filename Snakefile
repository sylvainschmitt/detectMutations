## Sylvain SCHMITT
## 21/07/2021

configfile: "config/config.dag.yml"

# libraries, = glob_wildcards(config["libdir"] + "/{library}_1.fq.gz")

rule all:
    input:
        ## reference ##
        # expand("results/reference/{reference}.fa", reference=config["reference"]), # ref
        ## reads ##
        # expand("results/{library}/{library}_{strand}.trimmed.paired.fq.gz", library=libraries, strand=["1", "2"]),
        ## alignments ##
        expand("results/alns/{library}.md.cram", library=config["leaf"]), # alns
        ## mutations ##
        # "results/leaf_nontrunk_mutations.sql",
        # "results/trunk_raw_mutations.sql",
        ## qc ##
        "results/multiqc_report.html"

# Rules #

## Reference ##
include: "rules/cp_reference.smk"
include: "rules/bwa_index.smk"

## Reads ##
# include: "rules/trimmomatic.smk"
# include: "rules/fastqc.smk"
include: "rules/multiqc.smk"

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

## Mutations ##
# include: "rules/strelka2.smk"
# include: "rules/bedtools_subtract.smk"
# include: "rules/strelka2tsv_leaf.smk"
# include: "rules/strelka2sql_leaf.smk"
# include: "rules/strelka2tsv_trunk.smk"
# include: "rules/strelka2sql_trunk.smk"
