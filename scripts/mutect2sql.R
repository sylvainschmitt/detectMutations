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
                      CHROM = col_character(),
                      REF = col_character(),
                      ALT = col_character(),
                      AS_SB_TABLE = col_logical(),
                      MBQ = col_character(),
                      MFRL = col_character(),
                      MMQ = col_character(),
                      STR = col_logical(),
                      PNOISE2 = col_logical(),
                      tumor = col_character(),
                      normal = col_character(),
                      caller = col_character()
                    )) %>% bind_rows()
vroom_write(mutations, path = csv, delim = ",")
rm(mutations)
csv_to_sqlite(csv_name = csv, db_name = sql, table_name = "mutations")
