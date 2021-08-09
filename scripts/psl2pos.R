tsvin <- snakemake@input[[1]]
psl <-  snakemake@input[[2]]
tsvout <-  snakemake@output[[1]]

# tsvin <- "data/bordeaux/3P_mutations.tsv"
# psl <-  "results/3P/mutations.psl"

library(tidyverse)

N <- 500

alns <- read_tsv(psl, skip = 5, 
                col_names = c("matches", "misMatches", "repMatches", "nCount", 
                              "qNumInsert", "qBaseInsert",
                              "tNumInsert", "tBaseInsert", "strand", 
                              "qName", "qSize", "qStart", "qEnd", 
                              "tName", "tSize", "tStart", "tEnd", 
                              "blockCount", "blockSizes", "qStarts", "tStarts"), 
                col_types = cols(blockSizes = col_character(), tStarts = col_character())) %>% 
  separate(qName, c("scf", "snvStart", "snvStop")) %>% 
  mutate(pos = as.numeric(snvStart) + 501) %>% 
  mutate(`Locus ID` = paste0(scf, "_", pos)) %>% 
  select(-snvStart, -snvStop, -scf, -pos)

read_tsv(tsvin) %>% 
  filter(`Mutation category` == "SM") %>% 
  left_join(alns) %>% 
  dplyr::select(-X2, -`Mutation category`, -`Chromosomal location`, -`f(alt)pool >0.5%`) %>% 
  dplyr::rename(locus = `Locus ID`, origin = `Origin of the mutation`) %>% 
  mutate(POS0 = 500+1) %>% # PSL 0-indexed !!
  filter(qStart <= POS0, qEnd >= POS0) %>% 
  group_by(Mutation) %>% 
  filter(matches == max(matches)) %>% 
  separate_rows(blockSizes, qStarts, tStarts, sep = ",", convert = T) %>% 
  filter(!is.na(qStarts)) %>% 
  filter(qStarts <= POS0) %>% 
  mutate(posSNV = tStarts + POS0 - qStarts) %>% 
  mutate(tEnd = tStarts + blockSizes) %>% 
  filter(posSNV <= tEnd) %>% 
  dplyr::select(locus, Mutation, origin, `f(alt)pool`, tName, posSNV) %>% 
  dplyr::rename(POS = posSNV, CHROM = tName) %>% 
  unique() %>% 
  write_tsv(tsvout)
