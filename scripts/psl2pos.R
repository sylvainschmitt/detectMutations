tsvin <- snakemake@input[[1]]
psl <-  snakemake@input[[2]]
tsvout <-  snakemake@output[[1]]

tsvin <- "results/mutations/napoleon_mutations.tsv"
psl <-  "results/mutations/napoleon_mutations.psl"

library(tidyverse)

snv0 <- read_tsv(tsvin) %>% 
  dplyr::rename(tumor = Mutation, Mutation = SNV, REF = Ref, ALT = Alt) %>% 
  dplyr::select(-Lower_genom, -Upper_genome) %>% 
  separate(Contig, c("Contig", "POS"), convert = T)
aln <- read_tsv(psl, skip = 5, 
                col_names = c("matches", "misMatches", "repMatches", "nCount", 
                              "qNumInsert", "qBaseInsert",
                              "tNumInsert", "tBaseInsert", "strand", 
                              "qName", "qSize", "qStart", "qEnd", 
                              "tName", "tSize", "tStart", "tEnd", 
                              "blockCount", "blockSizes", "qStarts", "tStarts"), 
                col_types = cols(blockSizes = col_character(), tStarts = col_character())) %>% 
  separate(qName, c("X1", "Contig", "X2", "snvStart", "snvStop")) %>% 
  select(-X1, -X2) %>% 
  mutate(Contig = gsub("1Quercusroburisolateech66", "", Contig)) %>% 
  mutate(POS = as.numeric(snvStart) + 501) %>% 
  select(-snvStart, -snvStop)

snv <- snv0 %>% 
  left_join(aln) %>% 
  dplyr::select(-Contig, -POS, -Context, -Position) %>% 
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
  dplyr::select(Mutation, tName, posSNV, tumor, REF, ALT) %>% 
  dplyr::rename(POS = posSNV, CHROM = tName) %>% 
  unique()

write_tsv(snv, tsvout)
