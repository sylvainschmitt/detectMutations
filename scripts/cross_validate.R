tsvs <- snakemake@input
out <-  snakemake@output[[1]]

library(tidyverse)

# tsvs <- c("results/mutations/Fagus_sylvatica_v3_mutations_on_Fagus_sylvatica_v3.tsv", 
#           "results/mutations/Fagus_sylvatica_v3_mutations_on_Fagus_sylvatica_mutant_v0.tsv", 
#           "results/mutations/Fagus_sylvatica_v3_mutations_on_Fagus_sylvatica_revertant_v0.tsv")

lapply(tsvs, function(tsv) 
         read_tsv(tsv) %>% mutate(tsv= tsv)) %>% 
  bind_rows() %>% 
  mutate(tsv = gsub("results/mutations/", "", tsv)) %>% 
  mutate(tsv = gsub(".tsv", "", tsv)) %>% 
  separate(tsv, c("mutations", "reference"), "_on_") %>% 
  mutate(value = 1) %>% 
  filter(reference != "Fagus_sylvatica_v3") %>% 
  reshape2::dcast(SNV + CHROM + POS + REF + ALT ~ reference) %>% 
  mutate_all(funs(ifelse(is.na(.), 0, .))) %>% 
  mutate(validation  = "all") %>% 
  mutate(validation = ifelse(Fagus_sylvatica_mutant_v0 == 0, "only revertant", validation)) %>% 
  mutate(validation = ifelse(Fagus_sylvatica_revertant_v0 == 0, "only mutant", validation)) %>% 
  dplyr::select(-Fagus_sylvatica_mutant_v0, -Fagus_sylvatica_revertant_v0) %>% 
  write_tsv(out)
