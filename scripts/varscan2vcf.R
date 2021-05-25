#input
varscanfile <- snakemake@input[[1]]
snpfile <- snakemake@input[[2]]
# output
vcffile <-  snakemake@output[[1]] 

# manual
# varscanfile <- "results/N100_R2_AF0.1_NR7000/varscan/N100_R2_AF0.1_NR7000.snp"
# snpfile <- file.path("results/reference/", yaml::read_yaml("config/config.dev.yml")$snps)
# vcffile <-  "results/N100_R2_AF0.1_NR7000/varscan/N100_R2_AF0.1_NR7000.vcf"

library(tidyverse)
library(Biostrings)

cat("##fileformat=VCFv4.1\n##fileDate=20210521\n##source=varscan\n#CHROM	POS	ID	REF	ALT\n",
    file = vcffile)

snps <-   read_tsv(snpfile, skip = 11) %>% 
  dplyr::rename(CHROM = `#CHROM`) %>% 
  select(CHROM, POS, ID, REF, ALT) %>% 
  mutate(ID = ".") %>% 
  as_tibble() %>% 
  mutate_all(as.character) %>% 
  mutate_at("POS", as.numeric)

read_tsv(varscanfile) %>% 
  mutate(ID = ".", CHROM	= chrom, POS = position, REF =ref,	ALT	= var) %>% 
  dplyr::select(CHROM,	POS,	ID,	REF,	ALT) %>% 
  anti_join(snps) %>% 
  write_tsv(file = vcffile, append = T, col_names = F)
  