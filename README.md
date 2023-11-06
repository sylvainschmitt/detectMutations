detect Mutations - Fruits
================
Sylvain Schmitt
September 29, 2023

- [Installation](#installation)
- [Usage](#usage)
  - [Locally](#locally)
  - [HPC](#hpc)
- [Workflow](#workflow)
  - [Reference](#reference)
  - [Reads](#reads)
  - [Alignments](#alignments)
  - [Variants](#variants)
- [Results](#results)

[`singularity` &
`snakemake`](https://github.com/sylvainschmitt/snakemake_singularity)
workflow to detect mutations in fruits from Angela and Sixto.

# Installation

- [x] Python ≥3.5
- [x] Snakemake ≥5.24.1
- [x] Golang ≥1.15.2
- [x] Singularity ≥3.7.3
- [x] This workflow

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
git checkout fruits
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

### [cp_reference](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/cp_reference.py)

- Tools: `cp`

### [bwa_index](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/bwa_index.py)

- Tools: [`BWA index`](http://bio-bwa.sourceforge.net/bwa.shtml)
- Singularity:
  oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest

### [samtools_faidx](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/samtools_faidx.py)

- Tools:
  [`samtools faidx`](http://www.htslib.org/doc/samtools-faidx.html)
- Singularity:
  oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [gatk_dict](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/gatk_dict.py)

- Tools:
  [`gatk CreateSequenceDictionary`](https://gatk.broadinstitute.org/hc/en-us/articles/360037422891-CreateSequenceDictionary-Picard-)
- Singularity: docker://broadinstitute/gatk

## Reads

*Trim reads.*

### [trimmomatic](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/trimmomatic.py)

- Tools:
  [`Trimmomatic`](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
- Singularity:
  oras://registry.forgemia.inra.fr/gafl/singularity/trimmomatic/trimmomatic:latest

## Alignments

*Align reads against reference, mark duplicated, and report alignment
quality.*

### [bwa_mem](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/bwa_mem.py)

- Tools: [`BWA mem`](http://bio-bwa.sourceforge.net/bwa.shtml)
- Singularity:
  oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest

### [samtools_view](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/samtools_view.py)

- Tools: [`Samtools view`](http://www.htslib.org/doc/samtools-view.html)
- Singularity:
  oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [samtools_sort](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/samtools_sort.py)

- Tools: [`Samtools sort`](http://www.htslib.org/doc/samtools-sort.html)
- Singularity:
  oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [samtools_index](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/samtools_index.py)

- Tools:
  [`Samtools index`](http://www.htslib.org/doc/samtools-index.html)
- Singularity:
  oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

## Variants

*Detect mutations.*

### [gatk_haplotypecaller](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/gatk_haplotypecaller.py)

- Tools:
  [`gatk HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
- Singularity: docker://broadinstitute/gatk

### [gatk_genomicdbimport](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/gatk_genomicdbimport.py)

- Tools:
  [`gatk GenomicsDBImport`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
- Singularity: docker://broadinstitute/gatk

### [gatk_genotypegvcfs](https://github.com/sylvainschmitt/detectMutations/blob/fruits/rules/ggatk_genotypegvcfs.py)

- Tools:
  [`gatk GenomicsDBImport`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
- Singularity: docker://broadinstitute/gatk

# Results
