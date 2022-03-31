callfile <- snakemake@input[[1]]
mutationsfile <- snakemake@input[[2]]
statsfile <-  snakemake@output[[1]] 

library(tidyverse)
library(vcfR)

mutations <- read_tsv(mutationsfile, col_types = cols(
  CHROM = col_character(),
  POS = col_double(),
  REF = col_character(),
  TYPE = col_character(),
  ALT = col_character()
)) %>% 
  mutate(True = 1) %>%
  filter(REF != "N")
call <- try(read.vcfR(callfile, verbose = F)@fix, silent = T)
if(!inherits(call, "try-error")){
  call <- as.data.frame(call) %>% 
    mutate_at(c("CHROM", "REF", "ALT"), as.character) %>% 
    mutate(POS = as.numeric(as.character(POS))) %>% 
    mutate(Called = 1) %>% 
    mutate(REF = substr(REF, 1, 1), ALT = substr(ALT, 1, 1))
} else {
  call <- data.frame(CHROM = character(), POS = double(), REF = character(), ALT = character(), Called = integer())
}
stats <- full_join(mutations, call, by = c("CHROM", "POS", "REF", "ALT")) %>% 
  mutate_at(c("True", "Called"), funs(ifelse(is.na(.), 0, .))) %>% 
  mutate(Confusion = recode(paste0(True, Called), "01" = "FP", "10" = "FN", "11" = "TP")) %>% 
  group_by(Confusion) %>% 
  summarise(N = n()) %>% 
  reshape2::dcast(mutationsfile + callfile ~ Confusion, value.var = "N")
if(!("FP" %in% names(stats)))
  stats$FP <- 0
if(!("FN" %in% names(stats)))
  stats$FN <- 0
if(!("TP" %in% names(stats)))
  stats$TP <- 0
stats <- mutate(stats, Precision = round(TP/(TP+FP), 2), Recall = round(TP/(TP+FN), 2)) %>% 
  mutate(callfile = callfile) %>% 
  mutate(mutationsfile = mutationsfile)
write_tsv(stats, statsfile)
