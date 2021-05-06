detect Mutations
================
Sylvain Schmitt
April 20, 2021

  - [Installation](#installation)
  - [Usage](#usage)
      - [Get data](#get-data)
      - [Locally](#locally)
      - [HPC](#hpc)
  - [Workflow](#workflow)
      - [Index](#index)
      - [Reads QC](#reads-qc)
      - [Trimming](#trimming)
      - [Mapping](#mapping)
      - [Recalibration](#recalibration)
      - [Alignments QC](#alignments-qc)
      - [Germline variant calling](#germline-variant-calling)
      - [Somatic variant calling](#somatic-variant-calling)
      - [Calls QC](#calls-qc)
  - [Results](#results)

[`singularity` &
`snakemake`](https://github.com/sylvainschmitt/snakemake_singularity)
workflow to detect mutations with several alignment and mutation
detection tools.

![Direct acyclic graph.](dag/dag.svg)

# Installation

  - [x] Python ≥3.5
  - [x] Snakemake ≥5.24.1
  - [x] Golang ≥1.15.2
  - [x] Singularity ≥3.7.3
  - [x] This workflow

<!-- end list -->

``` bash
# Python
sudo apt-get install python3.5
# Snakemake
sudo apt install snakemake`
# Golang
export VERSION=1.15.8 OS=linux ARCH=amd64  # change this as you need
wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz && \
sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
source ~/.bashrc
# Singularity
mkdir -p ${GOPATH}/src/github.com/sylabs && \
  cd ${GOPATH}/src/github.com/sylabs && \
  git clone https://github.com/sylabs/singularity.git && \
  cd singularity
git checkout v3.7.3
cd ${GOPATH}/src/github.com/sylabs/singularity && \
  ./mconfig && \
  cd ./builddir && \
  make && \
  sudo make install
# detect Mutations
git clone git@github.com:sylvainschmitt/detectMutations.git
cd detectMutations
```

# Usage

## Get data

*Generate data using the [generate
Mutations](https://github.com/sylvainschmitt/generateMutations)
workflow.*

``` bash
git clone git@github.com:sylvainschmitt/generateMutations.git
cd ../generateMutations
snakemake --use-singularity --cores 4
cd ../detectMutations
bash scripts/get_data.sh
```

## Locally

``` bash
snakemake -np # dry run
snakemake --dag | dot -Tsvg > dag/dag.svg # dag
snakemake --use-singularity --cores 4 # run
snakemake --use-singularity --cores 1 --verbose # debug
snakemake --report report.html # report
```

## HPC

``` bash
module purge ; module load bioinfo/snakemake-5.25.0 # for test on node
snakemake -np # dry run
sbatch job.sh ; watch 'squeue -u sschmitt' # run
less detMut.*.err # snakemake outputs, use MAJ+F
less detMut.*.out # snakemake outputs, use MAJ+F
snakemake --dag | dot -Tsvg > dag/dag.svg # dag
module purge ; module load bioinfo/snakemake-5.8.1 ; module load system/Python-3.6.3 # for report
snakemake --report report.html # report
```

# Workflow

## Index

*Index reference for software to work with.*

### [mv\_data](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/mv_data.smk)

  - Tools: `mv`

### [bwa\_index](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/bwa_index.smk)

  - Tools: [`BWA index`](http://bio-bwa.sourceforge.net/bwa.shtml)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest

### [samtools\_faidx](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/samtools_faidx.smk)

  - Tools: [`samtools
    faidx`](http://www.htslib.org/doc/samtools-faidx.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [gatk\_dict](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/gatk_dict.smk)

  - Tools: [`gatk
    CreateSequenceDictionary`](https://gatk.broadinstitute.org/hc/en-us/articles/360036729911-CreateSequenceDictionary-Picard-)
  - Singularity: docker://broadinstitute/gatk

## Reads QC

*Report read quality.*

### [fastqc](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/fastqc.smk)

  - Tools:
    [`fastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)
  - Singularity: docker://biocontainers/fastqc:v0.11.9\_cv8

### [multiqc](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/multiqc.smk)

  - Tools: [`MultiQC`](https://multiqc.info/)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/multiqc/multiqc:latest

## Trimming

*Trimming of low quality reads or bases is not implemented yet
([`Trimmomatic`](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)).*

## Mapping

*Align reads against reference.*

### [bwa\_mem](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/bwa_mem.smk)

  - Tools: [`BWA mem`](http://bio-bwa.sourceforge.net/bwa.shtml)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest

### [samtools\_sort](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/samtools_sort.smk)

  - Tools: [`Samtools
    sort`](http://www.htslib.org/doc/samtools-sort.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [samtools\_index](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/samtools_index.smk)

  - Tools: [`Samtools
    index`](http://www.htslib.org/doc/samtools-index.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [gatk\_markduplicates](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/gatk_markduplicates.smk)

  - Tools: [`gatk
    MarkDuplicates`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
  - Singularity: docker://broadinstitute/gatk

## Recalibration

*Recalibration necessitate confidence known sites. Without reference
SNPs DB, the user need to build one. Consequently this step is currently
skipped.*

## Alignments QC

*Report alignment quality.*

### [samtools\_stats](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/samtools_stats.smk)

  - Tools: [`Samtools
    stats`](http://www.htslib.org/doc/samtools-stats.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [qualimap](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/qualimap.smk)

  - Tools:
    [`QualiMap`](http://qualimap.conesalab.org/doc_html/command_line.html)
  - Singularity: docker://pegi3s/qualimap

## Germline variant calling

*Detect SNPs in the germline using `freebayes` and `gatk` (`strelka`,
`manta`, `tiddit`, etc can be further implemented).*

### [freebayes\_somatic](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/freebayes_somatic.smk)

  - Tools: [`freebayes`](https://github.com/freebayes/freebayes)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/freebayes/freebayes:latest

### GATK

#### [gatk\_haplotypecaller](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/gatk_haplotypecaller.smk)

  - Tools: [`gatk
    HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
  - Singularity: docker://broadinstitute/gatk

#### [gatk\_genotypegvcfs](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/gatk_genotypegvcfs.smk)

  - Tools: [`gatk
    GenotypeGVCFs`](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)
  - Singularity: docker://broadinstitute/gatk

## Somatic variant calling

*Detect mutations between germline and mutated somatic tissue using
`freebayes` and `Mutect2` (`strelka2`, `manta`, `ascat`,
`control-freec`, `msisensor`, etc can be further implemented).*

### [freebayes\_somatic](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/freebayes_somatic.smk)

  - Tools: [`freebayes`](https://github.com/freebayes/freebayes)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/freebayes/freebayes:latest

### [gatk\_mutect2](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/gatk_mutect2.smk)

  - Tools: [`gatk
    Mutect2`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
  - Singularity: docker://broadinstitute/gatk

## Calls QC

*Report calls quality.*

### [bcftools\_stats](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/bcftools_stats.smk)

  - Tools:
    [`bcftools_stats`](http://samtools.github.io/bcftools/bcftools.html#stats)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest

### [vcftools\_stats](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/vcftools_stats.smk)

  - Tools: [`vcftools --TsTv-by-count --TsTv-by-qual
    --FILTER-summary`](http://vcftools.sourceforge.net/man_latest.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/vcftools/vcftools:latest

<!-- template: -->

<!-- ### [command_sub](https://github.com/sylvainschmitt/detectMutations/blob/main/rules/command_sub.smk) -->

<!-- * Tools: [`command sub`](link) -->

<!-- * Singularity: link -->

# Results

> Note that `freebayes` has mistaken two SNPs at low coverage (\~75X).

<table>

<caption>

Generated mutations and their detection with different callers.

</caption>

<thead>

<tr>

<th style="text-align:left;">

Chromosome

</th>

<th style="text-align:right;">

Position

</th>

<th style="text-align:left;">

Reference

</th>

<th style="text-align:left;">

Alternative

</th>

<th style="text-align:left;">

Type

</th>

<th style="text-align:right;">

[`freebayes`](https://github.com/freebayes/freebayes)

</th>

<th style="text-align:right;">

[`Mutect2`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Qrob\_Chr01:0-1000

</td>

<td style="text-align:right;">

334

</td>

<td style="text-align:left;">

A

</td>

<td style="text-align:left;">

C

</td>

<td style="text-align:left;">

transversion1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr01:0-1000

</td>

<td style="text-align:right;">

342

</td>

<td style="text-align:left;">

C

</td>

<td style="text-align:left;">

G

</td>

<td style="text-align:left;">

transversion2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr01:0-1000

</td>

<td style="text-align:right;">

538

</td>

<td style="text-align:left;">

T

</td>

<td style="text-align:left;">

A

</td>

<td style="text-align:left;">

transversion2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr01:0-1000

</td>

<td style="text-align:right;">

609

</td>

<td style="text-align:left;">

A

</td>

<td style="text-align:left;">

C

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr01:0-1000

</td>

<td style="text-align:right;">

710

</td>

<td style="text-align:left;">

C

</td>

<td style="text-align:left;">

T

</td>

<td style="text-align:left;">

transition

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr01:0-1000

</td>

<td style="text-align:right;">

817

</td>

<td style="text-align:left;">

T

</td>

<td style="text-align:left;">

A

</td>

<td style="text-align:left;">

transversion2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr02:0-1000

</td>

<td style="text-align:right;">

352

</td>

<td style="text-align:left;">

G

</td>

<td style="text-align:left;">

C

</td>

<td style="text-align:left;">

transversion2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr02:0-1000

</td>

<td style="text-align:right;">

663

</td>

<td style="text-align:left;">

A

</td>

<td style="text-align:left;">

C

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr02:0-1000

</td>

<td style="text-align:right;">

665

</td>

<td style="text-align:left;">

G

</td>

<td style="text-align:left;">

C

</td>

<td style="text-align:left;">

transversion2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr02:0-1000

</td>

<td style="text-align:right;">

702

</td>

<td style="text-align:left;">

T

</td>

<td style="text-align:left;">

C

</td>

<td style="text-align:left;">

transition

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr02:0-1000

</td>

<td style="text-align:right;">

837

</td>

<td style="text-align:left;">

G

</td>

<td style="text-align:left;">

T

</td>

<td style="text-align:left;">

transversion1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Qrob\_Chr02:0-1000

</td>

<td style="text-align:right;">

984

</td>

<td style="text-align:left;">

G

</td>

<td style="text-align:left;">

A

</td>

<td style="text-align:left;">

transition

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

</tbody>

</table>
