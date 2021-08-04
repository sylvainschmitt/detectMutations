# input
vcf <- snakemake@input[[1]]
tsv <-  snakemake@output[[1]] 
tumor_sample <- snakemake@wildcards$tumor
normal_sample <- snakemake@wildcards$normal

# manual
# vcf <- "results/mutations/L1_vs_L2_mutect2.vcf"
# tsv <- "results/mutations/L1_vs_L2_mutect2.tsv"
# tumor_sample <- "L1"
# normal_sample <- "L2"


library(tidyverse)
library(vcfR)

raw <- read.vcfR(vcf, verbose = F) %>% 
  vcfR2tidy()

gt <- select(raw$gt, -gt_AF, -gt_F1R2, -gt_F2R1, -gt_GQ, -gt_GT, -gt_PGT, -gt_PID, -gt_PL, -gt_PS, -gt_SB, -gt_GT_alleles) %>% 
  separate(gt_AD, c("gt_refCount", "gt_altCount"), convert = T) %>% 
  arrange(ChromKey, POS) %>% 
  mutate(tissue = ifelse(Indiv == tumor_sample, "tumor", "normal")) %>% 
  select(-Indiv) %>% 
  mutate(gt_AF = gt_altCount/gt_DP) %>% 
  reshape2::melt(c("ChromKey", "POS", "tissue")) %>% 
  rowwise() %>% 
  mutate(variable = gsub("gt", tissue, variable)) %>% 
  ungroup() %>% 
  select(-tissue) %>% 
  reshape2::dcast(ChromKey + POS ~ variable) %>% 
  as.tbl()

select(raw$fix, -ID, -QUAL, -FILTER, -AS_UNIQ_ALT_READ_COUNT, -CONTQ, -GERMQ, -NCount, -OCM, -PON, -ROQ, -RPA, -RU, -SEQQ, -STRANDQ, -STRQ) %>% 
  filter(!grepl(",", ALT)) %>% 
  mutate_at(c("MPOS", "NALOD", "NLOD", "POPAF", "TLOD"), as.numeric) %>% 
  left_join(gt) %>% 
  select(-ChromKey) %>% 
  mutate(tumor = tumor_sample) %>% 
  mutate(normal = normal_sample) %>% 
  mutate(caller = "strelka2") %>% 
  write_tsv(tsv)
