detect Mutations
================
Sylvain Schmitt
April 20, 2021

  - [Align](#align)
      - [bwa](#bwa)
          - [bwa map](#bwa-map)
          - [PICARD SortSam](#picard-sortsam)
          - [samtools index](#samtools-index)
          - [samtools view](#samtools-view)
          - [samtools index](#samtools-index-1)
          - [samtools sort](#samtools-sort)
      - [GEM](#gem)
      - [Novoalign](#novoalign)
  - [Detect](#detect)
      - [Mutect2](#mutect2)
  - [Miscellaneous](#miscellaneous)
      - [Commands](#commands)
      - [Direct Acyclic Graph](#direct-acyclic-graph)
      - [Resources](#resources)

Development of a [`singularity` &
`snakemake`](https://github.com/sylvainschmitt/snakemake_singularity)
workflow to detect mutations with several alignment and mutation
detection tools.

# Align

*Align reads against reference.*

## bwa

### bwa map

*Align in sam.*

  - Rules:
    [`bwa_map`](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/bwa_map.smk)
  - Tools: [`BWA mem`](http://bio-bwa.sourceforge.net/bwa.shtml)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest

### PICARD SortSam

### samtools index

### samtools view

*Convert sam to bam.*

  - Rules:
    [`samtools_view`](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/samtools_view.smk)
  - Tools: [`samtools
    view`](http://www.htslib.org/doc/samtools-view.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### samtools index

### samtools sort

## GEM

## Novoalign

# Detect

*Detect and filter mutations in alignments.*

## Mutect2

  - Rules:
    [`mutect2`](https://github.com/sylvainschmitt/generateMutations/blob/main/rules/mutect2.smk)
  - Tools:
    [`Mutect2`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
  - Singularity: docker://alexcoppe/gatk

# Miscellaneous

## Commands

*To run locally.*

``` bash
snakemake -np 
snakemake --dag | dot -Tsvg > dag/dag.svg
snakemake --use-singularity --cores 4
snakemake --report report.html
```

*To run on HPC.*

``` bash
module purge ; module load bioinfo/snakemake-5.8.1 # for test on node
snakemake -np # dry run
sbatch job.sh ; watch 'squeue -u sschmitt' # run
less genMut.*.err # snakemake outputs, use MAJ+F
less genMut.*.out # snakemake outputs, use MAJ+F
snakemake --dag | dot -Tsvg > dag/dag.svg # dag
module purge ; module load bioinfo/snakemake-5.8.1 ; module load system/Python-3.6.3 # for report
snakemake --report report.html # report
```

## Direct Acyclic Graph

*Represent rules.*

![](dag/dag.svg)<!-- -->

## Resources

  - [TreeMutation
    pages](https://treemutation.netlify.app/mutations-detection.html#in-silico-mutations)
