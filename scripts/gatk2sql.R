tsv <- snakemake@input
csv <-  snakemake@output[[1]] 
sql <-  snakemake@output[[2]] 

# print(tsv)
# stop("test")

library(tidyverse)
library(vroom)
library(csv2sql)

mutations <- lapply(tsv, vroom, 
                    col_types = cols(
                      .default = col_double(),
                      caller = col_character(),	
                      CHROM = col_character(),	
                      REF = col_character(),	
                      ALT = col_character(),	
                      normal = col_character(),	
                      tumor = col_character()
                    )) %>% bind_rows()
vroom_write(mutations, path = csv, delim = ",")
rm(mutations)
csv_to_sqlite(csv_name = csv, db_name = sql, table_name = "mutations")
