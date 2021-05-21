#input
varscanfile <- snakemake@input[[1]]
# output
vcffile <-  snakemake@output[[1]] 

# manual
varscanfile <- "results/N100_R2_AF0.6_NR7000/varscan/N100_R2_AF0.6_NR7000.snp"
vcffile <-  "results/N100_R2_AF0.6_NR7000/varscan/N100_R2_AF0.6_NR7000.vcf"

library(tidyverse)
library(Biostrings)

cat("##fileformat=VCFv4.1\n##fileDate=20210521\n##source=varscan\n#CHROM	POS	ID	REF	ALT",
    file = vcffile)

read_tsv(varscanfile) %>% 
  mutate(ID = ".", CHROM	= chrom, POS = position, REF =ref,	ALT	= var) %>% 
  dplyr::select(CHROM,	POS,	ID,	REF,	ALT) %>% 
  write_tsv(file = vcffile, append = T, col_names = F)

##fileformat=VCFv4.1
##fileDate=20210521
##source=strelka
#CHROM	POS	ID	REF	ALT