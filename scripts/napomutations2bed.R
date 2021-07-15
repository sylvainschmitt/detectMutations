tsv <- snakemake@input[[1]] # "results/napoleon/napoleon_mutations.tsv"
fai <-  snakemake@input[[2]] # "results/napoleon/napoleon.fa.fai"
bed <-  snakemake@output[[1]] # "results/napoleon/napoleon_mutations.bed"

library(tidyverse)

N <- 500
mutations <- read_tsv(tsv) %>% 
  select(SNV, Contig) %>% 
  separate(Contig, c("Contig", "Pos"), convert = T) %>% 
  mutate(Start = Pos - N-1, Stop = Pos + N)
genome <- read_tsv(fai, 
                   col_names = c("Seq", "Size", "X1", "X2", "X3")) %>% 
  select(Seq) %>% 
  separate(Seq, paste0("X", 1:3), remove = F) %>% 
  mutate(Contig = gsub("1Quercusroburisolateech66", "", X2)) %>% 
  select(Seq, Contig)
genome %>% 
  left_join(mutations) %>% 
  filter(!is.na(SNV)) %>% 
  select(Seq, Start, Stop) %>% 
  write_tsv(bed, col_names = F)
