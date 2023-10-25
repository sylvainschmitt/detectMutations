configfile: "config/config.yml"

rule all:
    input:
        "results/angela.vcf",
        "results/sixto.vcf"
        # "results/angela_filtered.vcf",
        # "results/sixto_filtered.vcf"

# Rules #

## Reference & Reads ##
include: "rules/cp_reference.py"
include: "rules/bwa_index.py"
include: "rules/samtools_faidx.py"
include: "rules/gatk_dict.py"
include: "rules/trimmomatic.py"

## Alignments ##
include: "rules/bwa_mem.py"
include: "rules/samtools_view.py"
include: "rules/samtools_sort.py"
include: "rules/samtools_index.py"

## Variants ##
include: "rules/gatk_haplotypecaller.py"
include: "rules/gatk_genomicdbimport_angela.py"
include: "rules/gatk_genomicdbimport_sixto.py"
include: "rules/gatk_genotypegvcfs.py"
# include: "rules/bcftools_biallelic.py"
# include: "rules/gatk_snps.py"
