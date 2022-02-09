# snakemake
tsv <- snakemake@input[[1]]
fai <-  snakemake@input[[2]]
bed <-  snakemake@output[[1]]
N <- snakemake@params[['N']]
reference <- snakemake@params[['ref']]

library(tidyverse)

mutations <- read_tsv(tsv)
mutations <- mutations %>% 
  filter(reference == reference) %>% 
  mutate(SNV = paste0(CHROM, "#", POS)) %>% 
  select(SNV, CHROM, POS) %>% 
  mutate(Start = POS - N-1, Stop = POS + N)
genome <- read_tsv(fai, col_names = c("CHROM", "Size", "X1", "X2", "X3")) %>% 
  select(CHROM, Size)
left_join(genome, mutations) %>% 
  filter(!is.na(SNV)) %>% 
  mutate(Start = ifelse(Start < 0, 0, Start)) %>%
  mutate(Stop = ifelse(Stop > Size, Size, Stop)) %>%
  mutate(Start = ifelse(Start >= Stop, Stop - 1, Start)) %>% 
  select(CHROM, Start, Stop) %>% 
  write_tsv(bed, col_names = F)
