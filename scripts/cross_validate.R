tsvs <- snakemake@input
out <-  snakemake@output[[1]]

library(tidyverse)

lapply(tsvs, function(tsv) 
         read_tsv(tsv) %>% mutate(tsv= tsv)) %>% 
  bind_rows() %>% 
  mutate(tsv = gsub("results/mutations/", "", tsv)) %>% 
  mutate(tsv = gsub(".tsv", "", tsv)) %>% 
  separate(tsv, c("mutations", "reference"), "_on_") %>% 
  write_tsv(out)
