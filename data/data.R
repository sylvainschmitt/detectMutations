library(tidyverse)
read_tsv("filereport_read_run_PRJEB8388_tsv.txt") %>% 
  select(fastq_ftp) %>% 
  separate_rows(fastq_ftp, sep = ";") %>% 
  mutate(fastq_ftp = paste0("ftp://", fastq_ftp)) %>% 
  write_tsv("Bordeaux.list")
read_tsv("filereport_read_run_PRJNA327502_tsv.txt") %>% 
  select(fastq_ftp) %>% 
  separate_rows(fastq_ftp, sep = ";") %>% 
  mutate(fastq_ftp = paste0("ftp://", fastq_ftp)) %>% 
  write_tsv("Swiss.list")