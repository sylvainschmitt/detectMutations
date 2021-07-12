# input
raw_file <- snakemake@input[[1]]
# wildcards
normalAC <-  as.numeric(snakemake@config[["normalAC"]])
mutationsAC <- as.numeric(snakemake@config[["mutationsAC"]])
singleDP <- as.numeric(snakemake@config[["singleDP"]])
# output
filtered_file <-  snakemake@output[[1]] 

# # manual
# raw_file <- "results/mutations_unique/lower_vs_upper_on_Qrob_Chr01_strelka2.vcf"
# normalAC <- 0
# mutationsAC <- 10
# singleDP <- 100
# filtered_file <- "results/mutations_filtered/lower_vs_upper_on_Qrob_Chr01_strelka2.tsv"

# libraries
library(tidyverse)
library(vcfR)

# raw data
vcf <- read.vcfR(raw_file, verbose = F) %>% 
  vcfR2tidy()

# preparing gt
gt <- vcf$gt %>% 
  select(-gt_GT_alleles) %>% 
  reshape2::melt(c("ChromKey", "POS", "Indiv", "gt_DP", "gt_FDP", "gt_SDP", "gt_SUBDP"),
                 variable.name = "base", value.name = "count") %>% 
  mutate(base = gsub("gt_", "", base)) %>% 
  mutate(base = gsub("U", "", base)) %>% 
  left_join(select(vcf$fix, ChromKey, POS, REF, ALT)) %>% 
  rowwise() %>% 
  filter(base %in% c(REF, ALT)) %>% 
  mutate(base = ifelse(base == REF, "gt_refCounts", "gt_altCounts")) %>% 
  select(-REF, -ALT) %>% 
  reshape2::dcast(ChromKey+POS+Indiv+gt_DP+gt_FDP+gt_SDP+gt_SUBDP ~ base, value.var = "count") %>% 
  separate(gt_refCounts, c("gt_refCountT1", "gt_refCountT2"), convert = T) %>% 
  separate(gt_altCounts, c("gt_altCountT1", "gt_altCountT2"), convert = T) %>% 
  mutate(gt_AF = gt_altCountT1 / (gt_altCountT1 + gt_refCountT1))

# joining mutations infor
mut <- filter(gt, Indiv == "TUMOR") %>% select(-Indiv)
names(mut) <- gsub("gt", "mutation", names(mut))
norm <- filter(gt, Indiv == "NORMAL") %>% select(-Indiv)
names(norm) <- gsub("gt", "normal", names(norm))
mutations <- vcf$fix %>% 
  left_join(mut) %>% 
  left_join(norm)
rm(mut, norm)

# filtering
filtered <- mutations %>% 
  filter(normal_altCountT1 <= normalAC,
         normal_altCountT2 <= normalAC) %>% 
  filter(mutation_altCountT1 > mutationsAC,
         mutation_altCountT2 > mutationsAC) %>% 
  filter(normal_DP <= (singleDP*2), mutation_DP <= (singleDP*2), 
         normal_DP >= (singleDP/2), mutation_DP >= (singleDP/2)) %>% 
  mutate(tumor = snakemake@wildcards$tumor) %>% 
  mutate(normal = snakemake@wildcards$normal) %>% 
  mutate(caller = snakemake@wildcards$caller)
# to be added: selection of outputs

write_tsv(filtered, filtered_file)
