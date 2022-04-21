tsvin <- snakemake@input[[1]]
psl <-  snakemake@input[[2]]
tsvout <-  snakemake@output[[1]]

library(tidyverse)
library(vroom)

psl <- "results/mutations/Fagus_sylvatica_v3_mutations_on_Fagus_sylvatica_mutant_v0.psl"
tsvin <- "results/mutations/Fagus_sylvatica_v3_mutations.tsv"

aln <- vroom(psl, skip = 5, 
                    col_names = c("matches", "misMatches", "repMatches", "nCount", 
                                  "qNumInsert", "qBaseInsert",
                                  "tNumInsert", "tBaseInsert", "strand", 
                                  "qName", "qSize", "qStart", "qEnd", 
                                  "tName", "tSize", "tStart", "tEnd", 
                                  "blockCount", "blockSizes", "qStarts", "tStarts"), 
                    col_types = cols(blockSizes = col_character(), tStarts = col_character()))

aln <- sample_n(aln, 100)

aln <- aln %>% 
  separate(qName, c("CHROM", "StartStop"), ":") %>% 
  separate(StartStop, c("Start", "Stop"), "-", convert = T) %>% 
  left_join(read_tsv(tsvin)) %>% 
  mutate(POS = POS) %>% 
  mutate(POS0 = POS-Start) %>%
  filter(qStart <= POS0, qEnd >= POS0) %>% 
  group_by(SNV) %>% 
  filter(matches == max(matches)) %>% 
  separate_rows(blockSizes, qStarts, tStarts, sep = ",", convert = T) %>% 
  filter(!is.na(qStarts)) %>% 
  filter(qStarts <= POS0) %>% 
  mutate(posSNV = tStarts + POS0 - qStarts) %>% 
  mutate(tEnd = tStarts + blockSizes) %>% 
  filter(posSNV <= tEnd) %>% 
  dplyr::select(SNV, tName, posSNV, tumor, REF, ALT) %>% 
  dplyr::rename(POS = posSNV, CHROM = tName) %>% 
  unique()

vroom_write(aln, tsvout)
