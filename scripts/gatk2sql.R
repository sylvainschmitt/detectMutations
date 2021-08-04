tsv <- snakemake@input
csv <-  snakemake@output[[1]] 
sql <-  snakemake@output[[2]] 

# print(tsv)
# stop("test")
tsv <- list("results/mutations/upper_vs_lower_gatk.tsv",
            "results/mutations/lower_vs_upper_gatk.tsv")
csv <- "results/gatk_raw.csv"
sql <- "results/gatk_raw.sql"

library(tidyverse)
library(vroom)
library(csv2sql)


cols <- c("CHROM", "POS",
          "AC", "AF", "BaseQRankSum", "DP", "FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR",
          "gt_AD", "gt_DP", "gt_QG")

names(tsv) <- tsv
chromosomes <- vroom(tsv[1], col_names = cols) %>% 
  dplyr::select(CHROM) %>% 
  unique() %>% 
  unlist()

for(chromosome in chromosomes[1:12]){
  print(chromosome)
  lapply(tsv, vroom, 
         col_names = c("CHROM", "POS",
                       "AC", "AF", "BaseQRankSum", "DP", "FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR",
                       "gt_AD", "gt_DP", "gt_QG"),
         na = "."
  ) %>% 
    lapply(filter, CHROM == chromosome) %>% 
    bind_rows(.id = "sample") %>% 
    mutate(sample = gsub("results/mutations/", "", sample)) %>% 
    mutate(sample = gsub("_gatk.tsv", "", sample)) %>% 
    separate(sample, c("sample", "X2")) %>% 
    dplyr::select(-X2) %>% 
    separate(gt_AD, c("gt_refCount", "gt_altCount"), convert = T) %>% 
    arrange(CHROM, POS, desc(gt_altCount)) %>% 
    group_by(CHROM, POS) %>% 
    mutate(type = c("tumor", "normal")[1:n()]) %>% 
    mutate(tumor = sample[1]) %>% 
    mutate(normal = ifelse(sample[1] == "lower", "upper", "lower")) %>% 
    reshape2::melt(c("CHROM", "POS", "AC", "AF", "BaseQRankSum", "DP", "FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", "type", "sample", "tumor", "normal")) %>% 
    mutate(variable = gsub("gt", "", variable)) %>% 
    mutate(variable = paste0(type, variable)) %>% 
    dplyr::select(-type, -sample) %>%
    reshape2::dcast(CHROM + POS + tumor + normal ~ variable) %>%
    mutate(caller = "gatk") %>%
    arrange(CHROM, POS) %>% 
    mutate(normal_altCount = ifelse(is.na(normal_altCount), 0, normal_altCount)) %>% 
    mutate(tumor_AF = tumor_altCount/(tumor_altCount+tumor_refCount)) %>% 
    as.tbl() %>% 
    vroom_write(path = file.path("results", "tsvs", paste0(chromosome, ".csv")), delim = ",")
  gc()
}

print("unanchored")
lapply(tsv, vroom, 
       col_names = c("CHROM", "POS",
                     "AC", "AF", "BaseQRankSum", "DP", "FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR",
                     "gt_AD", "gt_DP", "gt_QG"),
       na = "."
) %>% 
  lapply(filter, !(CHROM %in% chromosomes[1:12])) %>% 
  bind_rows(.id = "sample") %>% 
  mutate(sample = gsub("results/mutations/", "", sample)) %>% 
  mutate(sample = gsub("_gatk.tsv", "", sample)) %>% 
  separate(sample, c("sample", "X2")) %>% 
  dplyr::select(-X2) %>% 
  separate(gt_AD, c("gt_refCount", "gt_altCount"), convert = T) %>% 
  arrange(CHROM, POS, desc(gt_altCount)) %>% 
  group_by(CHROM, POS) %>% 
  mutate(type = c("tumor", "normal")[1:n()]) %>% 
  mutate(tumor = sample[1]) %>% 
  mutate(normal = ifelse(sample[1] == "lower", "upper", "lower")) %>% 
  reshape2::melt(c("CHROM", "POS", "AC", "AF", "BaseQRankSum", "DP", "FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", "type", "sample", "tumor", "normal")) %>% 
  mutate(variable = gsub("gt", "", variable)) %>% 
  mutate(variable = paste0(type, variable)) %>% 
  dplyr::select(-type, -sample) %>%
  reshape2::dcast(CHROM + POS + tumor + normal ~ variable) %>%
  mutate(caller = "gatk") %>%
  arrange(CHROM, POS) %>% 
  mutate(normal_altCount = ifelse(is.na(normal_altCount), 0, normal_altCount)) %>% 
  mutate(tumor_AF = tumor_altCount/(tumor_altCount+tumor_refCount)) %>% 
  as.tbl() %>% 
  vroom_write(path = file.path("results", "tsvs", "unanchored.csv"), delim = ",")
gc()

lapply(list.files(file.path("results", "tsvs"), full.names = T), vroom) %>% 
  bind_rows() %>% 
  vroom_write(path = csv, delim = ",")

csv_to_sqlite(csv_name = csv, db_name = sql, table_name = "mutations")
