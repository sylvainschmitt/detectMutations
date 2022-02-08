# snakemake
raw_file <- snakemake@input[[1]]
base_file <- snakemake@output[[1]] 
robust_file <- snakemake@output[[2]] 

# libraries
library(tidyverse)
library(dbplyr)
library(dtplyr)
# DBI
# RSQLite

mutations <- DBI::dbConnect(RSQLite::SQLite(), raw_file) %>% tbl("mutations")
raw <- mutations %>% 
  collect() %>% 
  lazy_dt() %>%
  group_by(reference, CHROM, POS) %>% 
  filter(n() > 1) %>% 
  filter(normal_altCountT1 == min(normal_altCountT1)) %>% 
  filter(normal_altCountT1 == 0)
raw <- as_tibble(raw) %>% 
  mutate(reference = recode(reference,
                            "Fagus_sylvatica_v3"= "Ciron",
                            "Fagus_sylvatica_revertant_v0" = "Revertant" ,
                            "Fagus_sylvatica_mutant_v0" = "Mutant"
  ))
base <- filter(raw, 
               normal_DP <= 150, 
               normal_DP >= 25.5, 
               mutation_DP <= 150, 
               mutation_DP >= 25.5, 
               mutation_altCountT1 >= 5,
               mutation_AF < 0.5)
robust <- base[base$FILTER == "PASS",] 
write_tsv(base, file = base_file)
write_tsv(robust, file = robust_file)
