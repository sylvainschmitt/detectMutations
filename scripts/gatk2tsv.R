# input
vcf1 <- snakemake@input[[1]]
vcf2 <- snakemake@input[[2]]
tsv <-  snakemake@output[[1]] 

# # manual

# vcf1 <- "results/mutations/upper_vs_lower_on_Qrob_Chr01_gatk.vcf"
# vcf2 <- "results/mutations/lower_vs_upper_on_Qrob_Chr01_gatk.vcf"
# tsv <- "results/mutations/lower_vs_upper_on_Qrob_Chr01_gatk.tsv"

# libraries
library(tidyverse)
library(vcfR)
library(reshape2)

# raw data
cat("## Reading ...")
raw <- lapply(list(upper = vcf1, lower = vcf2), read.vcfR, verbose = F) %>% 
  lapply(vcfR2tidy)
meta <- raw$upper$meta
fix <- lapply(raw, `[[`, "fix")
gt <- lapply(raw, `[[`, "gt")
rm(raw)
cat(" done ##\n")

# gt
cat("## GT format ...")
gt <- gt %>% 
  bind_rows(.id = "branch") %>% 
  select(-ChromKey, -Indiv, -gt_GT, -gt_MIN_DP, -gt_PGT, -gt_PID, -gt_PL, 
         -gt_PS, -gt_RGQ, -gt_SB, -gt_GT_alleles) %>% 
  rename(AD = gt_AD, DP = gt_DP, GQ = gt_GQ) %>% 
  separate(AD, c("refCount", "altCount"), convert = T) %>% 
  group_by(POS) %>% 
  arrange(POS, altCount) %>% 
  mutate(type = NA) %>% 
  mutate(type = c("tumor", "normal")[1:n()])
cat(" done ##\n")  

# mutations
cat("## Mutations format ...")
mutations <- fix %>% 
  bind_rows(.id = "branch") %>% 
  select(-ChromKey, -ID, -FILTER, -END, -InbreedingCoeff, -RAW_MQandDP, -AC, -AF,
         -AN, -MLEAC, -MLEAF, -DP, -ExcessHet) %>% 
  left_join(gt, by = c("branch", "POS")) %>% 
  melt(c("type", "CHROM", "POS")) %>% 
  mutate(variable = paste0(variable, "_", type)) %>% 
  select(-type) %>% 
  reshape2::dcast(CHROM + POS ~ variable) %>% 
  select(-ALT_normal, -REF_normal, -QUAL_normal, -BaseQRankSum_normal,
         -FS_normal, -MQ_normal, -MQRankSum_normal, -QD_normal, 
         -ReadPosRankSum_normal, -SOR_normal) %>% 
  rename(ALT = ALT_tumor, REF = REF_tumor, QUAL = QUAL_tumor, BaseQRankSum = BaseQRankSum_tumor,
         FS = FS_tumor, MQ = MQ_tumor, MQRankSum = MQRankSum_tumor,
         QD = QD_tumor, ReadPosRankSum = ReadPosRankSum_tumor, SOR = SOR_tumor,
         normal = branch_normal, tumor = branch_tumor) %>% 
  mutate(caller = "gatk") %>% 
  select(caller, CHROM, POS, REF, ALT, normal, tumor, refCount_tumor, altCount_tumor, refCount_normal, altCount_normal,
         DP_tumor, DP_normal, GQ_tumor, GQ_normal, QUAL, BaseQRankSum, FS, MQ, MQRankSum, QD, ReadPosRankSum, SOR) %>% 
  mutate(normal = ifelse(is.na(normal) & tumor == "lower", "upper", normal)) %>% 
  mutate(normal = ifelse(is.na(normal) & tumor == "upper", "lower", normal)) %>% 
  mutate(altCount_normal = ifelse(is.na(altCount_normal), 0, altCount_normal))
cat(" done ##\n")  

cat("## Writing file ...")
write_tsv(mutations, tsv)
cat(" done ##\n")  
