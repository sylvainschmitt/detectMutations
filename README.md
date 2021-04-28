detect Mutations
================
Sylvain Schmitt
April 20, 2021

  - [Align](#align)
      - [Software 1](#software-1)
  - [Detect](#detect)
      - [Software 1](#software-1-1)
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

## Software 1

  - Rules:
    [`rule`](https://github.com/sylvainschmitt/generateMutations/blob/main/rules/rule.smk)
  - Tools: [`tool`](url)
  - Singularity: url
  - Parameter:

# Detect

*Detect and filter mutations in alignments.*

## Software 1

  - Rules:
    [`rule`](https://github.com/sylvainschmitt/generateMutations/blob/main/rules/rule.smk)
  - Tools: [`tool`](url)
  - Singularity: url
  - Parameter:

# Miscellaneous

## Commands

*To run locally.*

``` bash
snakemake -np 
snakemake --dag | dot -Tsvg > dag/dag.svg
snakemake --use-singularity --cores 4
snakemake --report results/report.html
```

## Direct Acyclic Graph

*Represent rules.*

## Resources

  - [TreeMutation
    pages](https://treemutation.netlify.app/mutations-detection.html#in-silico-mutations)
