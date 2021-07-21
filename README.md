detect Mutations - Angela
================
Sylvain Schmitt
Jully 21, 2021

  - [Sampling scheme](#sampling-scheme)
  - [Installation](#installation)
  - [Usage](#usage)
      - [Get data](#get-data)
      - [Locally](#locally)
      - [HPC](#hpc)
  - [Workflow](#workflow)
      - [Reference & reads](#reference--reads)
      - [Alignments](#alignments)
      - [Mutations](#mutations)
  - [Results](#results)

[`singularity` &
`snakemake`](https://github.com/sylvainschmitt/snakemake_singularity)
workflow to detect mutations corresponding to Angela sampling scheme.

![](dag/dag.minimal.svg)<!-- -->

# Sampling scheme

<img src="dag/sampling.png" width="636" />

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

Generate data using the [Angela branch of the generate
Mutations](https://github.com/sylvainschmitt/detectMutations/tree/angela)
workflow.

``` bash
git clone git@github.com:sylvainschmitt/generateMutations.git
git checkout -b angela
cd generateMutations
snakemake --use-singularity --cores 4
```

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

## Reference & reads

*Copy and index reference for software to work with and trim reads.*

### [samtools\_faidx\_split](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/samtools_faidx_split.smk)

  - Tools: [`samtools
    faidx`](http://www.htslib.org/doc/samtools-faidx.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest

### [bwa\_index](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/bwa_index.smk)

  - Tools: [`BWA index`](http://bio-bwa.sourceforge.net/bwa.shtml)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bwa/bwa:latest

### [trimmomatic](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/trimmomatic.smk)

  - Tools:
    [`Trimmomatic`](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/trimmomatic/trimmomatic:latest

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

## Mutations

### [strelka2](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/strelka2.smk)

  - Tools: [`Strelka2`](https://github.com/Illumina/strelka)
  - Singularity: docker://quay.io/wtsicgp/strelka2-manta

### [bedtools\_subtract](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/bedtools_subtract.smk)

  - Tools: [`bedtools
    subtract`](https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest

### [bedtools\_intersect](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/bedtools_intersect.smk)

  - Tools: [`bedtools
    intersect`](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)
  - Singularity:
    oras://registry.forgemia.inra.fr/gafl/singularity/bedtools/bedtools:latest

### [strelka2tsv](https://github.com/sylvainschmitt/detectMutations/blob/angela/rules/strelka2tsv.smk)

  - Script:
    [`strelka2tsv.R`](https://github.com/sylvainschmitt/detectMutations/blob/angela/scripts/strelka2tsv.R)
  - Singularity: to be added, currently uses local install

# Results

``` r
# tips
tips <- list.files("data/tips", pattern = ".tsv", full.names = T)
names(tips) <- gsub(".tsv", "", unlist(lapply(tips, basename)))
tips <- lapply(tips, read_tsv) %>% 
  bind_rows(.id = "sample") %>% 
  separate(sample, c("branch", "tip"), remove = F) %>% 
  mutate(origin = "tip")

# branches
branches <- list.files("data/branch", pattern = ".tsv", full.names = T)
names(branches) <- gsub(".tsv", "", unlist(lapply(branches, basename)))
branches <- lapply(branches, read_tsv) %>% 
  bind_rows(.id = "branch") %>% 
  rowwise %>% 
  mutate(tip = list(unique(tips$tip))) %>% 
  unnest() %>% 
  mutate(sample = paste0(branch, "_", tip)) %>% 
  mutate(origin = "branch")

# all
generated <- bind_rows(tips, branches) %>% 
  arrange(sample, origin) %>% 
  mutate(CHROM = gsub("_mutated", "", CHROM)) %>% 
  mutate(generated = 1)
rm(tips, branches)

# mutations
mutations <- list.files("results/mutations", pattern = ".tsv", full.names = T)
names(mutations) <- gsub("_on_Qrob_Chr01.tip.tsv", "", unlist(lapply(mutations, basename)))
mutations <- lapply(mutations, read_tsv) %>% 
  bind_rows(.id = "sample") %>% 
  unique() %>% 
  mutate(detected = 1)

# result
results <- generated %>% 
  full_join(mutations) %>% 
  mutate(generated = ifelse(is.na(generated), 0, generated)) %>% 
  mutate(detected = ifelse(is.na(detected), 0, detected)) %>% 
  mutate(status = recode(paste0(generated, detected), "11" = "TP", "10" = "FN", "01" = "FP"))
write_tsv(results, file = "experiment1.tsv")

# summary
results %>% 
  group_by(sample, status) %>% 
  summarise(N = n()) %>% 
  reshape2::dcast(sample ~ status) %>% 
  mutate(FN = 0, FP = 0) %>% 
  mutate(Precision = round(TP/(TP+FP), 2), 
         Recall = round(TP/(TP+FN), 2))
```
