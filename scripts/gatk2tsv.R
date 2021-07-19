# input
vcf <- snakemake@input[[1]]
tsv <-  snakemake@output[[1]] 

# # manual
# vcf <- "results/mutations_unique/lower_vs_upper_on_Qrob_Chr01_gatk.vcf"
# tsv <- "results/mutations_tsv/lower_vs_upper_on_Qrob_Chr01_gatk.tsv"

# loibraries
library(tidyverse)
library(vcfR)

# raw data
raw <- read.vcfR(vcf, verbose = F) %>% 
  vcfR2tidy()

# preparing mutations
mut <- raw$gt %>% 
  separate(gt_AD, c("gt_refCount", "gt_altCount"), convert = T)
names(mut) <- gsub("gt", "mutation", names(mut))
mutations <- raw$fix %>% 
  left_join(mut) %>% 
  mutate(tumor = snakemake@wildcards$tumor) %>% 
  mutate(normal = snakemake@wildcards$normal) %>% 
  mutate(caller = "gatk")
rm(mut)

write_tsv(mutations, tsv)
