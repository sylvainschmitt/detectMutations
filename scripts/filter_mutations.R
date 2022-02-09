infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]] 
lowDP <- snakemake@params[["lowDP"]] 
highDP <- snakemake@params[["highDP"]] 
minAC <- snakemake@params[["minAC"]] 
maxAF <- snakemake@params[["maxAFs"]] 

library(tidyverse)

DBI::dbConnect(RSQLite::SQLite(), infile) %>% 
  tbl("mutations") %>% 
  filter(
    normal_altCountT1 == 0,
    normal_DP <= highDP, 
    normal_DP >= lowDP, 
    mutation_DP <= highDP, 
    mutation_DP >= lowDP, 
    mutation_altCountT1 >= minAC,
    mutation_AF <= maxAF 
  ) %>% 
  collect() %>% 
  group_by(reference, CHROM, POS) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  mutate(Filter = ifelse(FILTER == "PASS", "robust", "base")) %>% 
  write_tsv(file = outfile)
