# snakemake
fai <- snakemake@input[[1]]
lib1 <- snakemake@input[[2]] 
lib2 <- snakemake@input[[3]] 
out <- snakemake@output[[1]] 
ref <- snakemake@params$reference
libs <- snakemake@params$library

# libraries
library(tidyverse)
library(circlize)
# vroom

# pars
n_chrom <- 15
depth <- 200

# data
genome <- read_tsv(fai, col_names = c("name", "length", "bytesindex", "basesperline", "bytesperline")) %>% 
  dplyr::select(name, length) %>% 
  mutate(start = 0, end = length) %>% 
  dplyr::select(-length) %>%  
  mutate(trsc = letters[1:nrow(.)], exon = LETTERS[1:nrow(.)]) %>% 
  arrange(desc(end)) %>% 
  top_n(n_chrom) %>% 
  as.data.frame
cov1 <- vroom::vroom(lib1, col_names = c("chrom", "start", "stop", "coverage")) %>% 
  filter(chrom %in% genome$name) %>% 
  mutate(coverage = ifelse(coverage > depth, depth, coverage)) %>% 
  as.data.frame()
cov2 <- vroom::vroom(lib2, col_names = c("chrom", "start", "stop", "coverage")) %>% 
  filter(chrom %in% genome$name) %>% 
  mutate(coverage = ifelse(coverage > depth, depth, coverage)) %>% 
  as.data.frame()

# circos
png(out, width = 1000, height = 1000)
circos.initializeWithIdeogram(genome)
circos.genomicTrack(cov1,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, col = "red", ...)
                    })
circos.genomicTrack(cov2,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, col = "blue", ...)
                    })
legend("bottomleft", pch = 16, col = c("blue", "red"), legend = libs, bty = "n")
text(0, 0, ref, cex = 1.5)
dev.off()
