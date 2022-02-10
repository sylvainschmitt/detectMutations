# snakemake
tsv <- snakemake@input[[1]]
fai <-  snakemake@input[[2]]
bed <-  snakemake@output[[1]]
tsvout <-  snakemake@output[[2]]
N <- snakemake@params[['N']]
ref <- snakemake@params[['ref']]
minAF <- snakemake@params[['minAF']]

library(tidyverse)

mutations <- read_tsv(tsv)
mutations <- mutations %>% 
  filter(reference == ref) %>% 
  filter(mutation_AF > minAF) %>% 
  filter(Filter == "robust") %>% 
  mutate(SNV = paste0(CHROM, "#", POS)) %>% 
  select(SNV, CHROM, POS, REF, ALT, tumor) %>% 
  mutate(Start = POS - N-1, Stop = POS + N)
genome <- read_tsv(fai, col_names = c("CHROM", "Size", "X1", "X2", "X3")) %>% 
  select(CHROM, Size)
mutations <- left_join(genome, mutations) %>% 
  filter(!is.na(SNV)) %>%
  rowwise() %>% 
  mutate(Start = max(0, Start)) %>% 
  mutate(Stop = min(Size, Stop))

write_tsv(select(mutations, CHROM, Start, Stop), bed, col_names = F)
write_tsv(select(mutations, SNV, CHROM, Start, POS, Stop, Size, REF, ALT, tumor), tsvout)
