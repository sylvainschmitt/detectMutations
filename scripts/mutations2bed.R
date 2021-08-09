tsv <- snakemake@input[[1]]
bed <-  snakemake@output[[1]]

# tsv <- "data/bordeaux/3P_mutations.tsv"

library(tidyverse)

N <- 500
read_tsv(tsv) %>% 
  filter(`Mutation category` == "SM") %>% 
  select(`Locus ID`) %>% 
  separate(`Locus ID`, c("scafold", "pos"), convert = T) %>% 
  mutate(start = pos - N-1, stop = pos + N) %>% 
  select(scafold, start, stop) %>% 
  write_tsv(bed, col_names = F)
