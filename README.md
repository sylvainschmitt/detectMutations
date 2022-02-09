detect Mutations - Angela
================
Sylvain Schmitt
Jully 21, 2021

  - [Installation](#installation)
  - [Usage](#usage)
      - [Locally](#locally)
      - [HPC](#hpc)
  - [Workflow](#workflow)
      - [Reference](#reference)
      - [Reads](#reads)
      - [Alignments](#alignments)
      - [Heterozygosity](#heterozygosity)
      - [Mutations](#mutations)

[`singularity` &
`snakemake`](https://github.com/sylvainschmitt/snakemake_singularity)
workflow to detect mutations corresponding to Angela sampling scheme.

![](dag/dag.minimal.svg)<!-- -->

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
git checkout angela
```

# Usage

## Locally

``` bash
snakemake -np -j 3 --resources mem_mb=10000 # dry run
snakemake --dag | dot -Tsvg > dag/dag.svg # dag
snakemake --use-singularity -j 3 --resources mem_mb=10000 # run
```

## HPC

``` bash
module load bioinfo/snakemake-5.25.0 # for test on node
snakemake -np # dry run
sbatch job.sh # run
snakemake --dag | dot -Tsvg > dag/dag.svg # dag
```

# Workflow

## Reference

*Copy and index reference.*

### [cp\_reference](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/cp_reference.smk)

  - Tools: `cp`

### [bwa\_index](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/bwa_index.smk)

  - Tools: [`BWA index`](http://bio-bwa.sourceforge.net/bwa.shtml)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest

## Reads

*Trim and quality check reads.*

### [trimmomatic](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/trimmomatic.smk)

  - Tools:
    [`Trimmomatic`](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/trimmomatic/trimmomatic:latest

### [fastqc](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/fastqc.smk)

  - Tools:
    [`fastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)
  - Singularity: docker://biocontainers/fastqc:v0.11.9\_cv8

### [multiqc](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/multiqc.smk)

  - Tools: [`MultiQC`](https://multiqc.info/)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/multiqc/multiqc:latest

## Alignments

*Align reads against reference, mark duplicated, and report alignment
quality.*

### [bwa\_mem](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/bwa_mem.smk)

  - Tools: [`BWA mem`](http://bio-bwa.sourceforge.net/bwa.shtml)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest

### [samtools\_view](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/samtools_view.smk)

  - Tools: [`Samtools
    view`](http://www.htslib.org/doc/samtools-view.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [samtools\_sort](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/samtools_sort.smk)

  - Tools: [`Samtools
    sort`](http://www.htslib.org/doc/samtools-sort.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [samtools\_index](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/samtools_index.smk)

  - Tools: [`Samtools
    index`](http://www.htslib.org/doc/samtools-index.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [gatk\_markduplicates](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/gatk_markduplicates.smk)

  - Tools: [`gatk
    MarkDuplicates`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
  - Singularity: docker://broadinstitute/gatk

### [samtools\_view\_md](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/samtools_view_md.smk)

  - Tools: [`Samtools
    view`](http://www.htslib.org/doc/samtools-view.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [samtools\_index\_md](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/samtools_index_md.smk)

  - Tools: [`Samtools
    index`](http://www.htslib.org/doc/samtools-index.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [samtools\_stats](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/samtools_stats.smk)

  - Tools: [`Samtools
    stats`](http://www.htslib.org/doc/samtools-stats.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [mosdepth](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/mosdepth.smk)

  - Tools: [`mosdepth`](https://github.com/brentp/mosdepth)
  - Singularity:
    docker://quay.io/biocontainers/mosdepth:0.2.4–he527e40\_0

## Heterozygosity

*Detect heterozygosity.*

### [gatk\_haplotypecaller](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/gatk_haplotypecaller.smk)

  - Tools: [`gatk
    HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
  - Singularity: docker://broadinstitute/gatk

### [gatk\_gathergvcfs](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/gatk_gathergvcfs.smk)

  - Tools: [`gatk
    GatherVcfs`](https://gatk.broadinstitute.org/hc/en-us/articles/360037422071-GatherVcfs-Picard-)
  - Singularity: docker://broadinstitute/gatk

### [gatk\_genomicsdbimport](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/gatk_genomicsdbimport.smk)

  - Tools: [`gatk
    GenomicsDBImport`](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport)
  - Singularity: docker://broadinstitute/gatk

### [gatk\_genotypegvcfs](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/gatk_genotypegvcfs.smk)

  - Tools: [`gatk
    GenotypeGVCFs`](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)
  - Singularity: docker://broadinstitute/gatk

### [gatk\_gathervcfs](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/gatk_gathervcfs.smk)

  - Tools: [`gatk
    GatherVcfs`](https://gatk.broadinstitute.org/hc/en-us/articles/360037422071-GatherVcfs-Picard-)
  - Singularity: docker://broadinstitute/gatk

### [bcftools\_biallelic](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/bcftools_biallelic.smk)

  - Tools:
    [`bcftools`](https://samtools.github.io/bcftools/bcftools.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest

### [gatk\_snps](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/gatk_snps.smk)

  - Tools: [`gatk
    VariantFiltration`](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)
  - Singularity: docker://broadinstitute/gatk

### [plink\_nonmissing](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/plink_nonmissing.smk)

  - Tools: [`plink`](https://www.cog-genomics.org/plink/)
  - Singularity:
    docker://quay.io/biocontainers/plink:1.90b6.21–h779adbc\_1

### [bcftools\_shared](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/bcftools_shared.smk)

  - Tools:
    [`bcftools`](https://samtools.github.io/bcftools/bcftools.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest

## Mutations

*Detect mutations in both cambiums and leaves.*

### [strelka2](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/strelka2.smk)

  - Tools: [`Strelka2`](https://github.com/Illumina/strelka)
  - Singularity: docker://quay.io/wtsicgp/strelka2-manta

### [bedtools\_subtract\_hz](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/bedtools_subtract_hz.smk)

  - Tools: [`bedtools
    subtract`](https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest

### [strelka2tsv](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/strelka2tsv.smk)

  - Script:
    [`strelka2tsv.R`](https://github.com/sylvainschmitt/detectMutations/blob/angela/scripts/strelka2tsv.R)
  - Singularity:
    <https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif>

### [strelka2tsv\_cambium](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/strelka2sql_cambium.smk)

  - Script:
    [`strelka2tsv.R`](https://github.com/sylvainschmitt/detectMutations/blob/angela/scripts/strelka2tsv.R)
  - Singularity:
    <https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif>

### [strelka2sql\_leaf](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/strelka2sql_leaf.smk)

  - Script:
    [`strelka2tsv.R`](https://github.com/sylvainschmitt/detectMutations/blob/angela/scripts/strelka2tsv.R)
  - Singularity:
    <https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif>

### [filter\_mutations](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/filter_mutations.smk)

  - Script:
    [`filter_mutations.R`](https://github.com/sylvainschmitt/detectMutations/blob/angela/scripts/filter_mutations.R)
  - Singularity:
    <https://github.com/sylvainschmitt/singularity-r-bioinfo/releases/download/0.0.3/sylvainschmitt-singularity-r-bioinfo.latest.sif>
