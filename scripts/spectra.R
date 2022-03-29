mutations_file <- snakemake@input[[1]]
genome_file <- snakemake@input[[2]]
genome_index_file <- snakemake@input[[3]]
out <-  snakemake@output[[1]]

library(tidyverse)
library(vroom)
library(Biostrings)
library(genomeIntervals)

# mutations_file <- "results/mutations/mutations.tsv"
# genome_file <- "results/reference/Fagus_sylvatica_mutant_v0.fa"
# genome_index_file <- "results/reference/Fagus_sylvatica_mutant_v0.fa.fai"

genome <- readDNAStringSet(genome_file)
genome_index <- read_tsv(genome_index_file, col_names = c("CHROM", "length", "bytesindex", "basesperline", "bytesperline")) %>% 
  dplyr::select(CHROM, length) %>% 
  mutate(start = 0, end = length) %>% 
  dplyr::select(-length)
names(genome) <- genome_index$CHROM
mutations <- vroom(mutations_file) %>% 
  filter(reference == snakemake@config[["reference"]])
seqs <- lapply(1:nrow(mutations), function(i)
  as.character(unlist(extractAt(genome[mutations$CHROM[i]], IRanges(mutations$POS[i]-1, width=3)))))
spectra <- mutations %>% 
  mutate(spectra = unlist(seqs)) %>% 
  mutate(REF5 = str_sub(spectra, 1, 1), REFspectra = str_sub(spectra, 2, 2), REF3 = str_sub(spectra, 3, 3)) %>% 
  dplyr::select(-REFspectra) %>% 
  mutate(type = paste0(REF, "->", ALT)) %>% 
  mutate(strand = recode(type,
                         "G->T" = "-",
                         "C->A" = "+",
                         "G->C" = "-",
                         "C->G" = "+",
                         "G->A" = "-",
                         "C->T" = "+",
                         "A->T" = "-",
                         "T->A" = "+",
                         "A->G" = "-",
                         "T->C" = "+",
                         "A->C" = "-",
                         "T->G" = "+"
  )) %>%
  mutate(type = recode(type,
                       "G->T" = "C->A",
                       "G->C" = "C->G",
                       "G->A" = "C->T",
                       "A->T" = "T->A",
                       "A->G" = "T->C",
                       "A->C" = "T->G"
  )) %>%
  mutate(class = ifelse(type %in% c("T->C", "C->T"), "transversion", "transition")) %>% 
  mutate(POS5 = REF5, POS3 = REF3) %>% 
  mutate(POS5r = recode(POS5, "A" = "T", "T" = "A", "C" = "G", "G" = "C")) %>%
  mutate(REFr = recode(REF, "A" = "T", "T" = "A", "C" = "G", "G" = "C")) %>%
  mutate(POS3r = recode(POS3, "A" = "T", "T" = "A", "C" = "G", "G" = "C")) %>%
  mutate(POS5f = ifelse(strand == "+", POS5, POS3r)) %>%
  mutate(REFf = ifelse(strand == "+", REF, REFr)) %>%
  mutate(POS3f = ifelse(strand == "+", POS3, POS5r)) %>% 
  mutate(spectra = paste0(POS5f, REFf, POS3f)) %>% 
  dplyr::select(CHROM, POS, REF, ALT, tumor, normal, mutation_AF, class, type, spectra)
write_tsv(spectra, out)

