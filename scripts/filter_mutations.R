# snakemake
infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]] 

# libraries
library(tidyverse)
library(dbplyr)
library(dtplyr)
# DBI
# RSQLite

DBI::dbConnect(RSQLite::SQLite(), infile) %>% 
  tbl("mutations") %>% 
  filter(
    tumor == "revertant",
    normal_altCountT1 == 0,
    normal_DP <= 150, 
    normal_DP >= 25.5, 
    mutation_DP <= 150, 
    mutation_DP >= 25.5, 
    mutation_altCountT1 >= 5,
    mutation_AF <= 0.5 
  ) %>% 
  collect() %>% 
  group_by(reference, CHROM, POS) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  mutate(Filter = ifelse(FILTER == "PASS", "robust", "base")) %>% 
  write_tsv(file = outfile)
