infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]] 
lowDP <- as.numeric(snakemake@params[["lowDP"]])
highDP <- as.numeric(snakemake@params[["highDP"]])
minAC <- as.numeric(snakemake@params[["minAC"]])
maxAF <- as.numeric(snakemake@params[["maxAF"]]) 

# infile <- "results/mutations_cambium/cambium_nonhz_mutations.sql"
# lowDP <- 50
# highDP <- 360
# minAC <- 10
# maxAF <- 0.5

library(tidyverse)

DBI::dbConnect(RSQLite::SQLite(), infile) %>% 
  tbl("mutations") %>% 
  collect() %>% 
  filter(
    normal_altCountT1 == 0,
    normal_DP <= highDP, 
    normal_DP >= lowDP, 
    mutation_DP <= highDP, 
    mutation_DP >= lowDP, 
    mutation_altCountT1 >= minAC,
    mutation_AF <= maxAF 
  ) %>% 
  group_by(tumor, CHROM, POS) %>% 
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(Filter = ifelse(FILTER == "PASS", "robust", "base")) %>%
  write_tsv(file = outfile)
